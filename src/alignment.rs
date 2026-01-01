use std::borrow::Cow;
use std::collections::BinaryHeap;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
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

// we implementes traits manually to avoid comparring the full record
// primary comparison is based using the random 'score'
// if there is tie (two read have the same score), then we use 'qname' to ensure the deterministic sorting
#[derive(Debug)]
struct ScoredRead {
    /// Read record from the alignment file
    record: rust_htslib::bam::Record,
    /// random score of the record
    score: u64,
}

impl ScoredRead {
    fn new(record: rust_htslib::bam::Record, score: u64) -> Self {
        Self { record, score }
    }
}

impl PartialEq for ScoredRead {
    fn eq(&self, other: &Self) -> bool {
        self.score == other.score && self.record.qname() == other.record.qname()
    }
}

impl Eq for ScoredRead {}

impl PartialOrd for ScoredRead {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

// if the random score is equal (very rare),
// then compare it based the qname record (which i think will be unique)
impl Ord for ScoredRead {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.score.cmp(&other.score) {
            std::cmp::Ordering::Equal => self.record.qname().cmp(other.record.qname()),
            std::cmp::Ordering::Greater => std::cmp::Ordering::Greater,
            std::cmp::Ordering::Less => std::cmp::Ordering::Less,
        }
    }
}

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

        // set up random number generator
        let mut rng = match self.seed {
            Some(s) => rand_pcg::Pcg64::seed_from_u64(s),
            None => {
                let seed = random();
                info!("Using seed: {seed}");
                rand_pcg::Pcg64::seed_from_u64(seed)
            }
        };

        // set up reader
        // because we use linear scan, i change IndexedReader to standard reader
        let mut reader =
            bam::Reader::from_path(&self.aln).context("Failed to read alignment file")?;

        let mut header = bam::Header::from_template(reader.header());

        // add rasusa program command line to header
        let program_record = self.program_entry(&header);
        header.push_record(&program_record);

        // set up header and writer
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

        // the sweep line algorithm starts here

        // the active set to store reads that currently in the scan
        // using max heap, it keeps the read with the highest score (worst priority) at the top.
        let target_depth = self.coverage as usize;
        let mut active_reads: BinaryHeap<ScoredRead> = BinaryHeap::with_capacity(target_depth);

        // we also pre alocated a vector here, and will be reused
        let mut survivors: Vec<ScoredRead> = Vec::with_capacity(target_depth);

        let mut current_tid = -1;
        let mut max_observed_depth: usize = 0;
        // track the last position so that we know whether the input is already sorted or no
        let mut last_pos: i64 = -1;

        // to get chromosome name for warning feature
        // only ask names if they actually exist
        let header_view = reader.header().clone();
        let chrom_names: Vec<String> = if header_view.target_count() > 0 {
            header_view
                .target_names()
                .iter()
                .map(|name_bytes| String::from_utf8_lossy(name_bytes).to_string())
                .collect()
        } else {
            Vec::new() // just give it an empty vector
        };

        // iterate through every read in the file
        for result in reader.records() {
            let record = result.context("Failed to parse BAM record")?;

            // skip unmapped reads (they do not have start and/or end coordinates)
            if record.is_unmapped() {
                continue;
            }

            // warning logic
            // tid ID maps to chromosome names, 0 means chrom 1
            let tid = record.tid();
            let pos = record.pos();

            // check if the input is not sorted yet
            if tid == current_tid && pos < last_pos {
                return Err(anyhow::anyhow!(
                    "Input is not sorted! Found read at pos {} after pos {} on chromosome {}",
                    pos,
                    last_pos,
                    tid
                ));
            }

            if tid != current_tid {
                // we just finished a chromosome. Did it ever reach target coverage?
                if current_tid != -1 && max_observed_depth < target_depth {
                    // Use .get() to safely check if the name exists.
                    // If chrom_names is empty (no_index.bam), this just skips the warning.
                    // but actually if no chromosome, means it is unmapped and will not go here?
                    if let Some(name) = chrom_names.get(current_tid as usize) {
                        warn!(
                            "Chromosome {} never reached the requested depth (Max observed: {}X)",
                            name, max_observed_depth
                        );
                    }
                } else if let Some(name) = chrom_names.get(current_tid as usize) {
                    info!("Subsampling completed successfully for: {name}");
                }
                // Reset for new chromosome
                current_tid = tid;
                max_observed_depth = 0;
                active_reads.clear();
            }

            let start = record.pos();

            // assign deterministic random priority
            let score: u64 = rng.gen();

            // clear and move everything from the heap to the preallocated vector
            // to check whether the reads in the heap already expire
            // the end pos of any read in the heap is lesser or equal than the current start pos.
            // i use drain method: https://doc.rust-lang.org/std/collections/struct.BinaryHeap.html#method.drain
            survivors.extend(active_reads.drain());

            for scored_read in survivors.drain(..) {
                if scored_read.record.reference_end() <= start {
                    // the reads survived!
                    // and we write it into the result
                    writer
                        .write(&scored_read.record)
                        .context("Failed to write record")?;
                } else {
                    // still active,
                    // put it back in heap
                    active_reads.push(scored_read);
                }
            }

            // note: survivors.drain(..) automatically clears the vector items but keeps the capacity.
            // with target depth usually 50 - 100 maybe, i think this operation is not very expensive.

            // when a new read begins:
            let new_read = ScoredRead::new(record, score);

            if active_reads.len() < target_depth {
                // if the heap has records fewer than N reads, just accept it.
                active_reads.push(new_read);
            } else {
                // if the heap is full, only accept it if its score is smaller than the current worst (highest score) in the heap
                let worse_score = active_reads.peek().unwrap().score;
                if new_read.score < worse_score {
                    active_reads.pop();
                    active_reads.push(new_read);
                }
            }

            // update the last position
            last_pos = pos;

            // track stats for Warning
            max_observed_depth = max_observed_depth.max(active_reads.len());
        }

        // warning logic: the last chromosome
        if current_tid != -1 && max_observed_depth < target_depth {
            // Use .get() to safely check if the name exists.
            // If chrom_names is empty (no_index.bam), this just skips the warning.
            if let Some(name) = chrom_names.get(current_tid as usize) {
                warn!(
                    "Chromosome {} never reached the requested depth (Max observed: {}X)",
                    name, max_observed_depth
                );
            }
        } else if let Some(name) = chrom_names.get(current_tid as usize) {
            info!("Subsampling completed successfully for: {name}");
        }

        // any remaining read in the heap survive
        // write into the subsampled result.
        for scored_read in active_reads {
            writer
                .write(&scored_read.record)
                .context("Failed to write record")?;
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
/// This function is from the original subsampling implementation, and kept for reference/testing
#[allow(dead_code)]
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
    use tempfile::NamedTempFile;

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

    // #[test]
    // fn bam_no_index_fails() {
    //     let infile = "tests/cases/no_index.bam";
    //     let passed_args = vec![SUB, infile, "-c", "1"];
    //     let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

    //     cmd.args(passed_args).assert().failure();
    // }

    // because it doesn't require index anymore
    #[test]
    fn bam_no_index_is_ok() {
        let infile = "tests/cases/no_index.bam";
        let passed_args = vec![SUB, infile, "-c", "1"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().success();
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

    #[test]
    fn test_scored_read_ordering() {
        let r1 = bam::Record::new();
        let small = ScoredRead {
            record: r1.clone(),
            score: 10,
        };

        let r2 = bam::Record::new();
        let big = ScoredRead {
            record: r2.clone(),
            score: 99,
        };

        assert!(big > small, "Higher score should be 'Greater' in ordering");
        assert!(small < big, "Lower score should be 'Less' in ordering");
    }

    #[test]
    fn test_scored_read_ordering_tie() {
        let mut r1 = bam::Record::new();
        r1.set_qname(b"readA");
        let record1 = ScoredRead {
            record: r1.clone(),
            score: 21,
        };

        let mut r2 = bam::Record::new();
        r2.set_qname(b"readB");
        let record2 = ScoredRead {
            record: r2.clone(),
            score: 21,
        };

        assert!(record2 > record1);
    }

    #[test]
    fn test_output_never_exceeds_target_depth() {
        use rust_htslib::bam::Read;

        let input_path = PathBuf::from("tests/cases/test.bam");
        let target_depth = 5;

        let output_file = NamedTempFile::new().unwrap();
        let output_path = output_file.path().to_path_buf();

        let mut align = Alignment {
            aln: input_path,
            output: Some(output_path.clone()),
            output_type: Some(bam::Format::Bam),
            coverage: target_depth,
            seed: Some(2109),
            step_size: 100, // we do not need this
        };

        // run the sweep line algorithm
        align.run().expect("Algorithm failed");

        // verify the resulted depth
        let mut reader = bam::Reader::from_path(&output_path).unwrap();

        let mut records: Vec<bam::Record> = reader.records().map(|r| r.unwrap()).collect();

        // sort by chromosome and then by position
        records.sort_by(|a, b| {
            let tid_cmp = a.tid().cmp(&b.tid());
            if tid_cmp == std::cmp::Ordering::Equal {
                a.pos().cmp(&b.pos())
            } else {
                tid_cmp
            }
        });

        let mut active_ends: Vec<i64> = Vec::new();
        let mut max_obs_depth = 0;
        let mut current_tid = -1;

        for record in records {
            let tid = record.tid();
            let start = record.pos();
            let end = record.reference_end();

            // new chromosome check
            if tid != current_tid {
                active_ends.clear();
                current_tid = tid;
            }

            // remove reads that ended before this one started, they are not contribute to the depth anymore
            active_ends.retain(|&e| e > start);

            // add current read
            active_ends.push(end);

            // check
            let current_depth = active_ends.len() as u32;
            if current_depth > max_obs_depth {
                max_obs_depth = current_depth;
            }
        }

        assert!(
            max_obs_depth <= target_depth,
            "The resulted depth exceeds the target depth"
        );
    }

    #[test]
    fn bam_is_not_sorted_fails() {
        let infile = "tests/cases/test_not_sorted.bam";
        let passed_args = vec![SUB, infile, "-c", "1"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().failure();
    }

    // helper function to run subsamping aln and get the resulted read names
    fn run_aln_get_reads_result(input: &Path, seed: Option<u64>) -> Vec<String> {
        let target_depth = 5;
        let out = NamedTempFile::new().unwrap();

        let mut aln1 = Alignment {
            aln: input.to_path_buf(),
            output: Some(out.path().to_path_buf()),
            output_type: Some(bam::Format::Bam),
            coverage: target_depth,
            seed,
            step_size: 100, // we do not need this anymore
        };

        aln1.run().expect("Subsampling failed");

        bam::Reader::from_path(out.path()).unwrap().records()
            .map(|r| {
                let record = r.unwrap();
                String::from_utf8_lossy(record.qname()).to_string()
            })
            .collect()
    }

    #[test]
    fn test_reproducibility_same_seed() {
        let input_path = Path::new("tests/cases/test.bam");
        let seed = Some(2109);

       let names1 = run_aln_get_reads_result(&input_path, seed);
       let names2 = run_aln_get_reads_result(&input_path, seed);

        // comapre the length and the read names
        assert_eq!(names1.len(), names2.len(), "Different read count");
        assert_eq!(names1, names2, "Different reads selected")
    }

    #[test]
    fn test_reproducibility_diff_seed() {
        let input_path = Path::new("tests/cases/test.bam");
        let seed1 = Some(21);
        let seed2 = Some(09);

       let names1 = run_aln_get_reads_result(&input_path, seed1);
       let names2 = run_aln_get_reads_result(&input_path, seed2);

        // comapre the length and the read names
        assert_ne!(names1, names2, "Same reads selected")
    }
}
