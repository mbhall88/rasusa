use std::borrow::Cow;
use std::collections::BinaryHeap;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;
use log::{info, warn};
use rand::{random, Rng, SeedableRng};
use rust_htslib::bam;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::{Format, Header, Read};

use crate::cli::check_path_exists;
use crate::Runner;

const RASUSA: &str = "rasusa";

// we implementes traits manually to avoid comparing the full record
// primary comparison is based using the random 'score'
// if there is tie (two read have the same score), then we use 'qname' to ensure the deterministic sorting still
#[derive(Debug)]
struct ScoredRead {
    /// alignment record (SAM/BAM/CRAM)
    record: rust_htslib::bam::Record,

    /// The deterministic score assigned to this record
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

    /// [DEPRECATED - Ignored in sweepline implementation]
    /// When a region has less than the desired coverage, the step size to move along the chromosome
    /// to find more reads.
    ///
    /// The lowest of the step and the minimum end coordinate of the reads in the region will be used.
    /// This parameter can have a significant impact on the runtime of the subsampling process.
    #[arg(long, default_value_t = 100, value_name = "INT", hide = true)] // hide from help
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

        // a helper closure for giving a warning or report after finished scanning a chromosome
        let depth_report = |tid: i32, max_depth: usize| {
            // only report if we processed a valid chromosome
            if tid == -1 {
                return;
            }
            // get the name for the chromsome, if it doesnt exist, give them unknown with askterisk (*)
            let name = chrom_names
                .get(tid as usize)
                .map(|n| n.as_str())
                .unwrap_or("*");

            // Did it ever reach target depth?
            if max_depth < target_depth {
                warn!(
                    "Chromosome {} never reached the requested depth (Max observed: {}X)",
                    name, max_depth
                );
            } else {
                info!("Subsampling complete for: {}", name);
            }
        };

        // iterate thriugh every read in the file
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
                depth_report(current_tid, max_observed_depth);
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

        // report a warning or info forn the last chromosome
        depth_report(current_tid, max_observed_depth);

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

        bam::Reader::from_path(out.path())
            .unwrap()
            .records()
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

        let names1 = run_aln_get_reads_result(input_path, seed);
        let names2 = run_aln_get_reads_result(input_path, seed);

        // comapre the length and the read names
        assert_eq!(names1.len(), names2.len(), "Different read count");
        assert_eq!(names1, names2, "Different reads selected")
    }

    #[test]
    fn test_reproducibility_diff_seed() {
        let input_path = Path::new("tests/cases/test.bam");
        let seed1 = Some(21);
        let seed2 = Some(9);

        let names1 = run_aln_get_reads_result(input_path, seed1);
        let names2 = run_aln_get_reads_result(input_path, seed2);

        // comapre the length and the read names
        assert_ne!(names1, names2, "Same reads selected")
    }

    #[test]
    fn unknown_input_extension_fails() {
        let input = Path::new("tests/cases/test");

        let mut aln = Alignment {
            aln: input.to_path_buf(),
            output: None,
            output_type: Some(bam::Format::Bam),
            coverage: 2,
            seed: Some(2109),
            step_size: 100, // we do not need this anymore
        };
        assert!(aln.run().is_err());
    }

    #[test]
    fn unknown_output_extension_fails() {
        let input = Path::new("tests/cases/test.bam");
        let output = Path::new("tests/cases/result");

        let mut aln = Alignment {
            aln: input.to_path_buf(),
            output: Some(output.to_path_buf()),
            output_type: None,
            coverage: 2,
            seed: Some(2109),
            step_size: 100, // we do not need this anymore
        };
        assert!(aln.run().is_err());
    }

    #[test]
    fn test_equality_record() {
        let mut r1 = bam::Record::new();
        r1.set_qname(b"readA");
        let record1 = ScoredRead {
            record: r1.clone(),
            score: 21,
        };

        let mut r2 = bam::Record::new();
        r2.set_qname(b"readA");
        let record2 = ScoredRead {
            record: r2.clone(),
            score: 21,
        };

        assert_eq!(record1, record2);
    }
}
