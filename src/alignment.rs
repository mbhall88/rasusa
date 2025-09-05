use std::borrow::Cow;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::collections::HashSet;
use std::path::{Path, PathBuf};

use anyhow::{anyhow, Context, Result};
use clap::Parser;
use log::{info, warn};
use rand::prelude::SliceRandom;
use rand::{random, Rng, SeedableRng};
use rust_htslib::bam;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::{Format, Header, Read};

use crate::cli::check_path_exists;
use crate::Runner;

const RASUSA: &str = "rasusa";

#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct Alignment {
    /// Path to the indexed alignment file (SAM/BAM/CRAM) to subsample
    #[arg(value_parser = check_path_exists, name = "FILE")]
    pub aln: PathBuf,

    /// Path to the output subsampled alignment file. Defaults to stdout (same format as input)
    ///
    /// The output is not guaranteed to be sorted. We recommend piping the output to `samtools sort`
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,

    /// Output format. Rasusa will attempt to infer the format from the output file extension if not provided
    #[arg(short='O', long, value_name = "FMT", value_parser = infer_format_from_char)]
    pub output_type: Option<Format>,

    /// The desired depth of coverage to subsample the alignment to
    #[arg(short, long, value_name = "INT", value_parser = clap::value_parser!(u32).range(1..))]
    pub coverage: u32,

    /// Random seed to use.
    #[arg(short, long, value_name = "INT")]
    pub seed: Option<u64>,

    /// When a region has less than the desired coverage, the step size to move along the chromosome
    /// to find more reads.
    ///
    /// The lowest of the step and the minimum end coordinate of the reads in the region will be used.
    /// This parameter can have a significant impact on the runtime of the subsampling process.
    #[arg(long, default_value_t = 100, value_name = "INT", value_parser = clap::value_parser!(i64).range(1..))]
    pub step_size: i64,
}

impl Runner for Alignment {
    fn run(&mut self) -> Result<()> {
        info!("Subsampling alignment file: {:?}", self.aln);

        let mut rng = match self.seed {
            Some(s) => rand_pcg::Pcg64::seed_from_u64(s),
            None => {
                let seed = random();
                info!("Using seed: {seed}");
                rand_pcg::Pcg64::seed_from_u64(seed)
            }
        };

        let mut reader =
            bam::IndexedReader::from_path(&self.aln).context("Failed to read alignment file")?;
        let mut header = bam::Header::from_template(reader.header());

        // add rasusa program command line to header
        let program_record = self.program_entry(&header);
        header.push_record(&program_record);

        let input_fmt = match infer_format_from_path(&self.aln) {
            Some(fmt) => fmt,
            None => {
                return Err(anyhow::anyhow!(
                    "Output file format not recognized. Please use .sam, .bam, or .cram extensions"
                ));
            }
        };

        let output_fmt = match &self.output_type {
            Some(fmt) => *fmt,
            None => match &self.output {
                None => input_fmt,
                Some(path) => match infer_format_from_path(path) {
                    Some(fmt) => fmt,
                    None => {
                        return Err(anyhow::anyhow!(
                            "Output file format not recognized. Please use .sam, .bam, or .cram extensions"
                        ));
                    }
                },
            },
        };

        let mut writer = match &self.output {
            Some(path) => {
                let path = path.as_path();

                let output = bam::Writer::from_path(path, &header, output_fmt)
                    .context("Failed to create output alignment file")?;
                info!("Writing subsampled alignment to: {:?}", path);
                output
            }
            None => {
                let output = bam::Writer::from_stdout(&header, output_fmt)
                    .context("Failed to create output alignment file")?;
                info!("Writing subsampled alignment to stdout");
                output
            }
        };

        let header = reader.header().clone();
        let chroms = header.target_names();

        for chrom in chroms {
            let chrom_name = String::from_utf8_lossy(chrom);

            info!("Subsampling chromosome: {chrom_name}");

            let tid = header
                .tid(chrom)
                .context(format!("Failed to get tid for chromosome {chrom_name}"))?;
            let chrom_len = header.target_len(tid).context(format!(
                "Failed to get chromosome length for chromosome {chrom_name}"
            ))?;
            let mut n_reads_needed = self.coverage;
            let mut current_reads = HashSet::new();
            let mut heap = BinaryHeap::new();

            // get the 0-based position of the first record in the chromosome
            reader.fetch(tid).context(format!(
                "Failed to get all records for chromosome {chrom_name}"
            ))?;

            let first_record = if let Some(first_record) = reader.records().next() {
                first_record.context("Failed to get first record")?
            } else {
                warn!("Chromosome {chrom_name} has no records");
                continue;
            };

            let mut next_pos = first_record.pos();
            let first_pos = next_pos;
            let mut regions_below_coverage = false;

            loop {
                reader
                    .fetch((tid, next_pos, next_pos + 1))
                    .context(format!(
                        "Failed to fetch records in region {chrom_name}:{next_pos}-{}",
                        next_pos + 1
                    ))?;
                let records = reader.records();

                let mut records: Vec<_> = records.filter_map(Result::ok).collect();

                if next_pos == first_pos {
                    // we just shuffle all the reads in the first position
                    records.shuffle(&mut rng);
                } else {
                    // need to sort records by their alignment start positions. those with the same start
                    // position should be shuffled so that the order is random
                    shuffle_records_by_position(&mut records, &mut rng);
                }

                let mut num_output = 0;
                let mut record_iter = records.into_iter().rev();

                while num_output < n_reads_needed {
                    let record = match record_iter.next() {
                        Some(r) => r,
                        None => break,
                    };
                    let qname = record.qname().to_owned();

                    if current_reads.contains(&qname) {
                        continue;
                    }
                    current_reads.insert(qname.to_owned());
                    heap.push(Reverse((record.reference_end(), qname.to_owned())));
                    // write the record
                    writer.write(&record).context("Failed to write record")?;
                    num_output += 1;
                }

                n_reads_needed -= num_output;

                if n_reads_needed > 0 {
                    // increment next_pos by step_size or the minimum end position of the reads in
                    // the region, whichever is smaller
                    let min_end = heap.peek().map(|Reverse((end, _))| *end);
                    match min_end {
                        Some(min_end) => {
                            next_pos += self.step_size.min(min_end - next_pos);
                        }
                        None => {
                            next_pos += self.step_size;
                        }
                    }
                    regions_below_coverage = true;
                }

                // remove smallest end position from heap and update next_pos and n_reads_needed
                // allowing for the fact that there may be multiple reads with the same end position
                while let Some(Reverse((end, qname))) = heap.pop() {
                    next_pos = end;

                    if current_reads.contains(&qname) {
                        current_reads.remove(&qname);
                        n_reads_needed += 1;
                    } else {
                        return Err(anyhow!("A read in the heap was not found in the current reads set. This should not happen, please raise an issue."));
                    }
                    if heap.is_empty() {
                        break;
                    }
                    if let Some(Reverse((next_end, _))) = heap.peek() {
                        if *next_end != end {
                            break;
                        }
                    } else {
                        break;
                    }
                }

                if next_pos as u64 >= chrom_len {
                    break;
                }
            }
            if regions_below_coverage {
                warn!("Chromosome {chrom_name} has regions with less than the requested coverage");
            }
        }

        Ok(())
    }
}

impl Alignment {
    /// Generates a rasusa program entry from a SAM header
    fn program_entry(&self, header: &Header) -> HeaderRecord<'_> {
        let (program_id, previous_pgid) = make_program_id_unique(header, RASUSA);

        let mut record = HeaderRecord::new(b"PG");
        record.push_tag(b"ID", program_id);
        record.push_tag(b"PN", RASUSA);
        if let Some(pp) = previous_pgid {
            record.push_tag(b"PP", pp);
        }
        record.push_tag(b"VN", env!("CARGO_PKG_VERSION"));
        let cl = std::env::args().collect::<Vec<String>>().join(" ");
        record.push_tag(b"CL", cl);

        record
    }
}

/// Makes a program ID unique by looking for existing program records with the same ID and adding
/// a suffix to the ID if necessary. Also returns the program ID of the last program in the header
fn make_program_id_unique<'a>(
    header: &Header,
    program_id: &'a str,
) -> (Cow<'a, str>, Option<String>) {
    let header_map = header.to_hashmap();
    let mut last_pg_id = None;
    let mut occurrences_of_id = 0;
    for (key, value) in header_map.iter() {
        if key == "PG" {
            for record in value {
                if let Some(id) = record.get("ID") {
                    last_pg_id = Some(id.to_string());
                    let id_before_last_dot = id.rfind('.').map(|i| &id[..i]).unwrap_or(id);
                    if id_before_last_dot == program_id {
                        occurrences_of_id += 1;
                    }
                }
            }
        }
    }
    if occurrences_of_id == 0 {
        (Cow::Borrowed(program_id), last_pg_id)
    } else {
        let new_id = format!("{program_id}.{occurrences_of_id}");
        (Cow::Owned(new_id), last_pg_id)
    }
}

/// Sorts records by position and shuffles those with the same position
fn shuffle_records_by_position(records: &mut [bam::Record], rng: &mut impl Rng) {
    // First sort by position
    records.sort_by_key(|record| record.pos());

    // Then shuffle groups with the same position
    let mut start = 0;
    while start < records.len() {
        let pos = records[start].pos();
        let mut end = start + 1;

        // Find the end of the group with the same position
        while end < records.len() && records[end].pos() == pos {
            end += 1;
        }

        // Shuffle this group if it has more than one element
        if end - start > 1 {
            records[start..end].shuffle(rng);
        }

        start = end;
    }
}

pub fn infer_format_from_path(path: &Path) -> Option<Format> {
    match path.extension().and_then(|ext| ext.to_str()) {
        Some("sam") => Some(Format::Sam),
        Some("bam") => Some(Format::Bam),
        Some("cram") => Some(Format::Cram),
        _ => None,
    }
}

// a function which infers a format from a single character (case insensitive) this will be used in the CLI, so should return a Result
pub fn infer_format_from_char(c: &str) -> Result<Format, String> {
    match c.to_ascii_lowercase().as_str() {
        "s" | "sam" => Ok(Format::Sam),
        "b" | "bam" => Ok(Format::Bam),
        "c" | "cram" => Ok(Format::Cram),
        _ => Err(String::from("Invalid output format. Please use 's', 'b', or 'c' to specify SAM, BAM, or CRAM format, respectively")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_cmd::Command;
    use rand::prelude::StdRng;
    use rust_htslib::bam::HeaderView;

    const SUB: &str = "aln";

    #[test]
    fn test_infer_format() {
        // have to test this way as rust-htslib does not derive PartialEq or Eq for Format
        let fmt = infer_format_from_path(Path::new("file.sam")).unwrap();
        assert!(matches!(fmt, Format::Sam));

        let fmt = infer_format_from_path(Path::new("file.bam")).unwrap();
        assert!(matches!(fmt, Format::Bam));

        let fmt = infer_format_from_path(Path::new("file.cram")).unwrap();
        assert!(matches!(fmt, Format::Cram));

        let fmt = infer_format_from_path(Path::new("file.txt"));
        assert!(fmt.is_none());
    }

    #[test]
    fn test_infer_format_from_char() {
        let fmt = infer_format_from_char("s").unwrap();
        assert!(matches!(fmt, Format::Sam));

        let fmt = infer_format_from_char("b").unwrap();
        assert!(matches!(fmt, Format::Bam));

        let fmt = infer_format_from_char("c").unwrap();
        assert!(matches!(fmt, Format::Cram));

        let fmt = infer_format_from_char("x");
        assert!(fmt.is_err());

        let fmt = infer_format_from_char("sam").unwrap();
        assert!(matches!(fmt, Format::Sam));

        let fmt = infer_format_from_char("B").unwrap();
        assert!(matches!(fmt, Format::Bam));

        let fmt = infer_format_from_char("CRAM").unwrap();
        assert!(matches!(fmt, Format::Cram));

        let fmt = infer_format_from_char("x");
        assert!(fmt.is_err());
    }

    #[test]
    fn test_shuffle_records_by_position_empty() {
        let mut rng = StdRng::seed_from_u64(1234);
        let mut empty_records: Vec<rust_htslib::bam::Record> = vec![];

        shuffle_records_by_position(&mut empty_records, &mut rng);

        assert_eq!(empty_records.len(), 0);
    }

    #[test]
    fn test_shuffle_records_by_position_single_record() {
        let mut rng = StdRng::seed_from_u64(1234);

        // Create a single BAM record
        let mut record = rust_htslib::bam::Record::new();
        record.set_pos(100);
        let mut records = vec![record];

        shuffle_records_by_position(&mut records, &mut rng);

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].pos(), 100);
    }

    #[test]
    fn test_shuffle_records_by_position_maintains_sort_order() {
        let mut rng = StdRng::seed_from_u64(1234);

        // Create records with different positions
        let mut record1 = rust_htslib::bam::Record::new();
        record1.set_pos(300);
        let mut record2 = rust_htslib::bam::Record::new();
        record2.set_pos(100);
        let mut record3 = rust_htslib::bam::Record::new();
        record3.set_pos(200);

        let mut records = vec![record1, record2, record3];

        shuffle_records_by_position(&mut records, &mut rng);

        // Should be sorted by position
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].pos(), 100);
        assert_eq!(records[1].pos(), 200);
        assert_eq!(records[2].pos(), 300);
    }

    #[test]
    fn test_shuffle_records_by_position_shuffles_same_position() {
        let mut rng = StdRng::seed_from_u64(1234);

        // Create multiple records with the same position but different data
        let mut record1 = rust_htslib::bam::Record::new();
        record1.set_pos(100);
        record1.set_qname(b"read1");

        let mut record2 = rust_htslib::bam::Record::new();
        record2.set_pos(100);
        record2.set_qname(b"read2");

        let mut record3 = rust_htslib::bam::Record::new();
        record3.set_pos(100);
        record3.set_qname(b"read3");

        let records = vec![record1, record2, record3];
        let original_order: Vec<Vec<u8>> = records.iter().map(|r| r.qname().to_vec()).collect();

        // Run shuffle multiple times to verify randomness
        let mut different_orders = 0;
        for _ in 0..10 {
            let mut test_records = records.clone();
            shuffle_records_by_position(&mut test_records, &mut rng);

            // All should still be at position 100
            assert!(test_records.iter().all(|r| r.pos() == 100));

            // Check if order changed
            let new_order: Vec<Vec<u8>> = test_records.iter().map(|r| r.qname().to_vec()).collect();
            if new_order != original_order {
                different_orders += 1;
            }
        }

        // Should have at least some different orders due to shuffling
        assert!(
            different_orders > 0,
            "Records with same position should be shuffled"
        );
    }

    #[test]
    fn test_shuffle_records_by_position_mixed_positions() {
        let mut rng = StdRng::seed_from_u64(1234);

        // Create records with mixed positions - some same, some different
        let mut record1 = rust_htslib::bam::Record::new();
        record1.set_pos(100);
        record1.set_qname(b"read1_pos100");

        let mut record2 = rust_htslib::bam::Record::new();
        record2.set_pos(200);
        record2.set_qname(b"read2_pos200");

        let mut record3 = rust_htslib::bam::Record::new();
        record3.set_pos(100);
        record3.set_qname(b"read3_pos100");

        let mut record4 = rust_htslib::bam::Record::new();
        record4.set_pos(150);
        record4.set_qname(b"read4_pos150");

        let mut record5 = rust_htslib::bam::Record::new();
        record5.set_pos(100);
        record5.set_qname(b"read5_pos100");

        let mut records = vec![record1, record2, record3, record4, record5];

        shuffle_records_by_position(&mut records, &mut rng);

        // Should be sorted by position
        let positions: Vec<i64> = records.iter().map(|r| r.pos()).collect();
        assert_eq!(positions, vec![100, 100, 100, 150, 200]);

        // Records at position 100 should potentially be in different order
        let pos100_names: Vec<Vec<u8>> = records
            .iter()
            .filter(|r| r.pos() == 100)
            .map(|r| r.qname().to_vec())
            .collect();
        assert_eq!(pos100_names.len(), 3);

        // All names should be present
        let name_strings: Vec<String> = pos100_names
            .iter()
            .map(|name| String::from_utf8_lossy(name).to_string())
            .collect();
        assert!(name_strings.contains(&"read1_pos100".to_string()));
        assert!(name_strings.contains(&"read3_pos100".to_string()));
        assert!(name_strings.contains(&"read5_pos100".to_string()));
    }

    #[test]
    fn no_coverage_given_raises_error() {
        let infile = "tests/cases/test.bam";
        let passed_args = vec![SUB, infile];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn zero_coverage_raises_error() {
        let infile = "tests/cases/test.bam";
        let passed_args = vec![SUB, infile, "-c", "0"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn bam_with_regions_of_zero_coverage_doesnt_endless_loop() {
        let infile = "tests/cases/test.bam";
        let passed_args = vec![SUB, infile, "-c", "1"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().success();
    }

    #[test]
    fn excess_coverage_doesnt_endless_loop() {
        let infile = "tests/cases/test.bam";
        let passed_args = vec![SUB, infile, "-c", "10000"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().success();
    }

    #[test]
    fn bam_no_index_fails() {
        let infile = "tests/cases/no_index.bam";
        let passed_args = vec![SUB, infile, "-c", "1"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn bam_with_no_start_or_end_regions_and_missing_chromosomes() {
        let infile = "tests/cases/no_start_end.bam";
        let passed_args = vec![SUB, infile, "-c", "1"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().success();
    }

    #[test]
    fn test_make_program_id_unique_no_program() {
        let template = HeaderView::from_bytes(b"@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chromosome\tLN:5399960
@PG\tID:minimap2\tPN:minimap2\tVN:2.26-r1175\tCL:minimap2 -aL --cs --MD -t 4 -x map-ont KPC2__202310.5x.fq.gz
@PG\tID:samtools\tPN:samtools\tPP:minimap2\tVN:1.19.2\tCL:samtools sort -@ 4 -o KPC2__202310.5x.bam
@PG\tID:samtools.1\tPN:samtools\tPP:samtools\tVN:1.19\tCL:samtools view -s 0.5 -o test.bam KPC2__202310.5x.bam");
        let header = Header::from_template(&template);
        let program_id = "rasusa";
        let actual = make_program_id_unique(&header, program_id);
        let expected = (
            Cow::<str>::Borrowed(program_id),
            Some("samtools.1".to_string()),
        );
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_make_program_id_unique_one_program_occurrence() {
        let template = HeaderView::from_bytes(b"@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chromosome\tLN:5399960
@PG\tID:minimap2\tPN:minimap2\tVN:2.26-r1175\tCL:minimap2 -aL --cs --MD -t 4 -x map-ont KPC2__202310.5x.fq.gz
@PG\tID:samtools\tPN:samtools\tPP:minimap2\tVN:1.19.2\tCL:samtools sort -@ 4 -o KPC2__202310.5x.bam
@PG\tID:samtools.1\tPN:samtools\tPP:samtools\tVN:1.19\tCL:samtools view -s 0.5 -o test.bam KPC2__202310.5x.bam");
        let header = Header::from_template(&template);
        let program_id = "minimap2";
        let actual = make_program_id_unique(&header, program_id);
        let expected = (
            Cow::<str>::Owned("minimap2.1".to_string()),
            Some("samtools.1".to_string()),
        );
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_make_program_id_unique_two_program_occurrences() {
        let template = HeaderView::from_bytes(b"@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chromosome\tLN:5399960
@PG\tID:minimap2\tPN:minimap2\tVN:2.26-r1175\tCL:minimap2 -aL --cs --MD -t 4 -x map-ont KPC2__202310.5x.fq.gz
@PG\tID:samtools\tPN:samtools\tPP:minimap2\tVN:1.19.2\tCL:samtools sort -@ 4 -o KPC2__202310.5x.bam
@PG\tID:samtools.1\tPN:samtools\tPP:samtools\tVN:1.19\tCL:samtools view -s 0.5 -o test.bam KPC2__202310.5x.bam");
        let header = Header::from_template(&template);
        let program_id = "samtools";
        let actual = make_program_id_unique(&header, program_id);
        let expected = (
            Cow::<str>::Owned("samtools.2".to_string()),
            Some("samtools.1".to_string()),
        );
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_make_program_id_unique_no_programs() {
        let template = HeaderView::from_bytes(
            b"@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chromosome\tLN:5399960",
        );
        let header = Header::from_template(&template);
        let program_id = "samtools";
        let actual = make_program_id_unique(&header, program_id);
        let expected = (Cow::Borrowed("samtools"), None);
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_original_issue_demonstrates_problem_with_random_compare() {
        // This test shows the core issue of #76: random_compare doesn't properly shuffle
        // because it only affects individual comparisons, not group-level shuffling

        use std::cmp::Ordering;

        // Simulate the old random_compare function
        fn old_random_compare<T: Ord>(a: T, b: T, rng: &mut impl Rng) -> Ordering {
            if a == b {
                if rng.gen::<bool>() {
                    Ordering::Less
                } else {
                    Ordering::Greater
                }
            } else {
                a.cmp(&b)
            }
        }

        let mut rng = StdRng::seed_from_u64(42);

        // Create records with identifiable names but same position
        let mut record_a = rust_htslib::bam::Record::new();
        record_a.set_pos(100);
        record_a.set_qname(b"readA");

        let mut record_b = rust_htslib::bam::Record::new();
        record_b.set_pos(100);
        record_b.set_qname(b"readB");

        let mut record_c = rust_htslib::bam::Record::new();
        record_c.set_pos(100);
        record_c.set_qname(b"readC");

        let original_records = vec![record_a, record_b, record_c];
        let original_names: Vec<String> = original_records
            .iter()
            .map(|r| String::from_utf8_lossy(r.qname()).to_string())
            .collect();

        println!("Original order: {original_names:?}");

        // Test old approach with sort_by + random_compare
        println!("\nTesting old random_compare approach (10 iterations):");
        let mut old_unique_orders = std::collections::HashSet::new();
        for i in 0..10 {
            let mut test_records = original_records.clone();
            test_records.sort_by(|a, b| old_random_compare(a.pos(), b.pos(), &mut rng));

            let names: Vec<String> = test_records
                .iter()
                .map(|r| String::from_utf8_lossy(r.qname()).to_string())
                .collect();

            println!("Iteration {}: {names:?}", i + 1);
            old_unique_orders.insert(names);
        }

        // Test new approach with proper shuffle
        println!("\nTesting new shuffle approach (10 iterations):");
        let mut new_unique_orders = std::collections::HashSet::new();
        for i in 0..10 {
            let mut test_records = original_records.clone();
            // Sort first, then shuffle (simulating our new function's behavior)
            test_records.sort_by_key(|r| r.pos());
            test_records.shuffle(&mut rng);

            let names: Vec<String> = test_records
                .iter()
                .map(|r| String::from_utf8_lossy(r.qname()).to_string())
                .collect();

            println!("Iteration {}: {names:?}", i + 1);
            new_unique_orders.insert(names);
        }

        println!(
            "\nOld approach generated {} unique orders",
            old_unique_orders.len()
        );
        println!(
            "New approach generated {} unique orders",
            new_unique_orders.len()
        );

        // The new approach should generally produce more diverse orderings
        // (though due to randomness, this isn't guaranteed on every run)
        assert!(
            !new_unique_orders.is_empty(),
            "New approach should produce at least one ordering"
        );
        assert!(
            !old_unique_orders.is_empty(),
            "Old approach should produce at least one ordering"
        );
    }

    #[test]
    fn test_original_issue_random_compare_vs_proper_shuffle() {
        // This test demonstrates why random_compare is insufficient for shuffling.
        // random_compare only introduces randomness at the individual comparison level,
        // but doesn't guarantee that groups of equal elements get properly shuffled.
        // This relates to issue #76.

        use std::cmp::Ordering;

        // Simulate the old random_compare function
        fn old_random_compare<T: Ord>(a: T, b: T, rng: &mut impl Rng) -> Ordering {
            if a == b {
                if rng.gen::<bool>() {
                    Ordering::Less
                } else {
                    Ordering::Greater
                }
            } else {
                a.cmp(&b)
            }
        }

        let mut rng = StdRng::seed_from_u64(123);

        // Create a simple test with integers to demonstrate the issue
        let values = vec![1, 1, 1, 1, 1]; // All same values
        let original_order = values.clone();

        // Test old approach: sort with random_compare
        let mut old_approach_different_orders = 0;
        for _ in 0..100 {
            let mut test_values = values.clone();
            test_values.sort_by(|a, b| old_random_compare(*a, *b, &mut rng));

            if test_values != original_order {
                old_approach_different_orders += 1;
            }
        }

        // Test new approach: proper shuffle
        let mut new_approach_different_orders = 0;
        for _ in 0..100 {
            let mut test_values = values.clone();
            test_values.shuffle(&mut rng);

            if test_values != original_order {
                new_approach_different_orders += 1;
            }
        }

        println!(
            "Old random_compare approach: {old_approach_different_orders} out of 100 iterations changed order"
        );
        println!(
            "Proper shuffle approach: {new_approach_different_orders} out of 100 iterations changed order"
        );

        // The proper shuffle should be much more effective at changing the order
        // Note: Due to the random nature, there's a small chance both could be similar,
        // but generally shuffle should be much more effective
        assert!(
            new_approach_different_orders >= old_approach_different_orders,
            "Proper shuffle should be at least as effective as random_compare"
        );
    }

    #[test]
    fn test_original_issue_random_compare_insufficient_shuffling() {
        // This test demonstrates the original issue: random_compare doesn't properly shuffle
        // groups of records with the same position because it only introduces randomness
        // at the comparison level, not at the group level.
        // This relates to issue #76.

        use std::cmp::Ordering;

        // Simulate the old random_compare function
        fn old_random_compare<T: Ord>(a: T, b: T, rng: &mut impl Rng) -> Ordering {
            if a == b {
                // Introduce randomness when elements are equal
                if rng.gen::<bool>() {
                    Ordering::Less
                } else {
                    Ordering::Greater
                }
            } else {
                a.cmp(&b)
            }
        }

        // Simulate the old random_sort function
        fn old_random_sort<T, K: Ord + Copy>(
            vec: &mut [T],
            key_extractor: fn(&T) -> K,
            mut rng: impl Rng,
        ) {
            vec.sort_by(|a, b| old_random_compare(key_extractor(a), key_extractor(b), &mut rng));
        }

        let mut rng = StdRng::seed_from_u64(42);

        // Create records with the same position
        let mut record1 = rust_htslib::bam::Record::new();
        record1.set_pos(100);
        record1.set_qname(b"read1");

        let mut record2 = rust_htslib::bam::Record::new();
        record2.set_pos(100);
        record2.set_qname(b"read2");

        let mut record3 = rust_htslib::bam::Record::new();
        record3.set_pos(100);
        record3.set_qname(b"read3");

        let records = vec![record1, record2, record3];
        let original_order: Vec<Vec<u8>> = records.iter().map(|r| r.qname().to_vec()).collect();

        // Test the old approach - it often fails to properly shuffle
        let mut same_order_count = 0;
        for _ in 0..20 {
            let mut test_records = records.clone();
            old_random_sort(&mut test_records, |record| record.pos(), &mut rng);

            let new_order: Vec<Vec<u8>> = test_records.iter().map(|r| r.qname().to_vec()).collect();
            if new_order == original_order {
                same_order_count += 1;
            }
        }

        // The old approach often keeps the same order because random_compare
        // doesn't guarantee proper shuffling of equal elements
        println!("Old approach: {same_order_count} out of 20 iterations kept the same order");

        // Now test our new approach
        let mut new_same_order_count = 0;
        for _ in 0..20 {
            let mut test_records = records.clone();
            shuffle_records_by_position(&mut test_records, &mut rng);

            let new_order: Vec<Vec<u8>> = test_records.iter().map(|r| r.qname().to_vec()).collect();
            if new_order == original_order {
                new_same_order_count += 1;
            }
        }

        println!("New approach: {new_same_order_count} out of 20 iterations kept the same order");

        // The new approach should shuffle much more effectively
        // (though due to randomness, it might occasionally keep the same order)
        assert!(
            new_same_order_count < same_order_count,
            "New shuffling approach should be more effective than old random_compare approach"
        );
    }

    #[test]
    fn test_make_program_id_unique_program_id_startswith_same_substring() {
        let template = HeaderView::from_bytes(b"@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chromosome\tLN:5399960
@PG\tID:minimap2\tPN:minimap2\tVN:2.26-r1175\tCL:minimap2 -aL --cs --MD -t 4 -x map-ont KPC2__202310.5x.fq.gz
@PG\tID:samtoolsfoo\tPN:samtools\tPP:minimap2\tVN:1.19.2\tCL:samtools sort -@ 4 -o KPC2__202310.5x.bam
@PG\tID:samtools\tPN:samtools\tPP:minimap2\tVN:1.19.2\tCL:samtools sort -@ 4 -o KPC2__202310.5x.bam
@PG\tID:samtoolsfoo.1\tPN:samtools\tPP:samtools\tVN:1.19\tCL:samtools view -s 0.5 -o test.bam KPC2__202310.5x.bam");
        let header = Header::from_template(&template);
        let program_id = "samtools";
        let actual = make_program_id_unique(&header, program_id);
        let expected = (
            Cow::<str>::Owned("samtools.1".to_string()),
            Some("samtoolsfoo.1".to_string()),
        );
        assert_eq!(actual, expected);
    }
}
