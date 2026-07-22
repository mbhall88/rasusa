mod args;
mod fetch;
mod header;
mod io;
mod model;
mod stream;
mod util;

pub use args::{Alignment, SubsamplingStrategy};
pub use header::{make_program_id_unique, program_entry};

use std::collections::HashSet;

use anyhow::{Context, Result};
use log::info;
use noodles::sam::Header;
use noodles_util::alignment;
use rustc_hash::FxBuildHasher;

use crate::Runner;
use io::AlignmentWriter;
use util::extract_name;

/// Set of read (QNAME) bytes, keyed with a fast non-cryptographic hasher since these are
/// populated/queried heavily in paired-end mate recovery and offer no adversarial-input risk.
pub(super) type NameSet = HashSet<Vec<u8>, FxBuildHasher>;

impl Runner for Alignment {
    fn run(&mut self) -> Result<()> {
        match self.strategy {
            SubsamplingStrategy::Stream => self.run_stream(),
            SubsamplingStrategy::Fetch => self.run_fetch(),
        }
    }
}

impl Alignment {
    // a function which infers data type (paired and or single end) by looking at the first 10 records
    fn check_pair(&self) -> Result<bool> {
        // set up reader
        let mut reader = alignment::io::reader::Builder::default()
            .build_from_path(&self.aln)
            .context("Failed to read alignment file")?;

        let header = reader.read_header()?;

        // read records only to check what type of the data
        for result in reader.records(&header).take(10) {
            let record = result.context("Failed to parse BAM record")?;
            if record.flags()?.is_segmented() {
                return Ok(true);
            }
        }
        Ok(false)
    }

    // a helper function to calculate target depth based on their input
    fn get_target_depth(&self, is_paired: bool) -> u32 {
        if is_paired {
            // divide the original target by 2 because we will add the mates back later
            (self.coverage / 2).max(1) // in case we hit 0
        } else {
            self.coverage
        }
    }

    // a helper function to scan the file linearly to find mates, only relevent on paired end illumina data
    fn recover_mates(
        &self,
        survivor_names: &mut NameSet,
        header: &Header,
        writer: &mut AlignmentWriter,
    ) -> Result<()> {
        info!("Recovering mates (last segment records)");
        let mut reader = alignment::io::reader::Builder::default()
            .build_from_path(&self.aln)
            .context("Failed to read alignment file")?;

        let _ = reader.read_header()?; // skip header

        for result in reader.records(header) {
            let record = result.context("Failed to parse BAM record")?;

            // only care about the last segment (read 2)
            if !record.flags()?.is_last_segment() {
                continue;
            }

            let qname: Vec<u8> = extract_name(&record);

            if survivor_names.contains(&qname) {
                writer
                    .write_record(header, &record)
                    .context("Failed to write mate records")?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_cmd::Command;
    use noodles_util::alignment::io::Format;
    use std::path::{Path, PathBuf};
    use tempfile::NamedTempFile;

    const SUB: &str = "aln";

    #[test]
    fn no_coverage_given_raises_error() {
        let infile = "tests/cases/test.bam";
        let passed_args = vec![SUB, infile];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn zero_coverage_raises_error_stream() {
        let infile = "tests/cases/test.bam";
        let passed_args = vec![SUB, infile, "-c", "0"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn zero_coverage_raises_error_fetch() {
        let infile = "tests/cases/test.bam";
        let passed_args = vec![
            SUB,
            infile,
            "-c",
            "0",
            "--strategy",
            "fetch",
            "--step-size",
            "5000",
        ];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn bam_with_regions_of_zero_coverage_doesnt_endless_loop_stream() {
        let infile = "tests/cases/test.bam";
        let passed_args = vec![SUB, infile, "-c", "1"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().success();
    }

    #[test]
    fn bam_with_regions_of_zero_coverage_doesnt_endless_loop_fetch() {
        let infile = "tests/cases/test.bam";
        let passed_args = vec![SUB, infile, "-c", "1", "--strategy", "fetch"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().success();
    }

    #[test]
    fn excess_coverage_doesnt_endless_loop_stream() {
        let infile = "tests/cases/test.bam";
        let passed_args = vec![SUB, infile, "-c", "10000"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().success();
    }

    #[test]
    fn excess_coverage_doesnt_endless_loop_fetch() {
        let infile = "tests/cases/test.bam";
        let passed_args = vec![SUB, infile, "-c", "10000", "--strategy", "fetch"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().success();
    }

    // because it doesn't require index anymore
    #[test]
    fn bam_no_index_is_ok_stream() {
        let infile = "tests/cases/no_index.bam";
        let passed_args = vec![SUB, infile, "-c", "1"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().success();
    }

    #[test]
    fn bam_no_index_fails_fetch() {
        let infile = "tests/cases/no_index.bam";
        let passed_args = vec![SUB, infile, "-c", "1", "--strategy", "fetch"];
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
    fn bam_is_not_sorted_fails() {
        let infile = "tests/cases/test_not_sorted.bam";
        let passed_args = vec![SUB, infile, "-c", "1"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().failure();
    }

    // helper function to run subsamping aln and get the resulted read names
    fn run_aln_get_reads_result(
        input: &Path,
        seed: Option<u64>,
        strategy: SubsamplingStrategy,
    ) -> Vec<String> {
        let target_depth = 3;
        let out = NamedTempFile::new().unwrap();

        let mut aln1 = Alignment {
            aln: input.to_path_buf(),
            output: Some(out.path().to_path_buf()),
            output_format: Some(Format::Bam),
            coverage: target_depth,
            seed,
            strategy,
            swap_distance: 5,
            step_size: 100,
            batch_size: 10_000,
        };

        aln1.run().expect("Subsampling failed");

        let mut reader = alignment::io::reader::Builder::default()
            .build_from_path(out.path())
            .unwrap();
        let header = reader.read_header().unwrap();

        reader
            .records(&header)
            .map(|r| {
                let rec = r.unwrap();
                // Get name as String from the record
                String::from_utf8_lossy(rec.name().unwrap()).to_string()
            })
            .collect()
    }

    #[test]
    fn test_reproducibility_stream_same_seed() {
        let input_path = Path::new("tests/cases/test.bam");
        let seed = Some(2109);

        let names1 = run_aln_get_reads_result(input_path, seed, SubsamplingStrategy::Stream);
        let names2 = run_aln_get_reads_result(input_path, seed, SubsamplingStrategy::Stream);

        // comapre the length and the read names
        assert_eq!(names1.len(), names2.len(), "Different read count");
        assert_eq!(names1, names2, "Different reads selected")
    }

    #[test]
    fn test_reproducibility_fetch_same_seed() {
        let input_path = Path::new("tests/cases/test.bam");
        let seed = Some(2109);

        let names1 = run_aln_get_reads_result(input_path, seed, SubsamplingStrategy::Fetch);
        let names2 = run_aln_get_reads_result(input_path, seed, SubsamplingStrategy::Fetch);

        // comapre the length and the read names
        assert_eq!(names1.len(), names2.len(), "Different read count");
        assert_eq!(names1, names2, "Different reads selected")
    }

    #[test]
    fn test_reproducibility_stream_diff_seed() {
        let input_path = Path::new("tests/cases/test.bam");
        let seed1 = Some(21);
        let seed2 = Some(9);

        let names1 = run_aln_get_reads_result(input_path, seed1, SubsamplingStrategy::Stream);
        let names2 = run_aln_get_reads_result(input_path, seed2, SubsamplingStrategy::Stream);

        // comapre the length and the read names
        assert_ne!(names1, names2, "Same reads selected")
    }

    #[test]
    fn test_reproducibility_fetch_diff_seed() {
        let input_path = Path::new("tests/cases/test.bam");
        let seed1 = Some(21);
        let seed2 = Some(9);

        let names1 = run_aln_get_reads_result(input_path, seed1, SubsamplingStrategy::Fetch);
        let names2 = run_aln_get_reads_result(input_path, seed2, SubsamplingStrategy::Fetch);

        // comapre the length and the read names
        assert_ne!(names1, names2, "Same reads selected")
    }

    #[test]
    fn unknown_input_extension_fails() {
        let input = Path::new("tests/cases/test");

        let mut aln = Alignment {
            aln: input.to_path_buf(),
            output: None,
            output_format: Some(Format::Bam),
            coverage: 2,
            seed: Some(2109),
            strategy: SubsamplingStrategy::Stream,
            swap_distance: 5,
            step_size: 100, // not used
            batch_size: 10_000,
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
            output_format: None,
            coverage: 2,
            seed: Some(2109),
            strategy: SubsamplingStrategy::Stream,
            swap_distance: 5,
            step_size: 100, // not used
            batch_size: 10_000,
        };
        assert!(aln.run().is_err());
    }

    #[test]
    fn test_paired_end_retention_stream() {
        let input_path = PathBuf::from("tests/cases/test.paired.bam");

        let target_depth: u32 = 2;
        let output = NamedTempFile::new().unwrap();

        // subsample to 2x depth
        let mut align = Alignment {
            aln: input_path,
            output: Some(output.path().to_path_buf()),
            output_format: Some(Format::Bam),
            coverage: target_depth,
            seed: Some(2109),
            strategy: SubsamplingStrategy::Stream,
            swap_distance: 5,
            step_size: 100,
            batch_size: 1000,
        };

        align.run().expect("Subsampling failed");

        // verify
        let mut reader = alignment::io::reader::Builder::default()
            .build_from_path(output.path())
            .unwrap();
        // read the header, get the length of chromosome
        let header = reader.read_header().unwrap();

        let mut r1_names: Vec<String> = Vec::new();
        let mut r2_names: Vec<String> = Vec::new();

        for result in reader.records(&header) {
            let record = result.unwrap();
            let name = record.name().unwrap().to_string();
            let flags = record.flags().unwrap();

            if flags.is_first_segment() {
                r1_names.push(name);
            } else if flags.is_last_segment() {
                r2_names.push(name);
            }
        }
        // make sure r1 and r2 names are identical because they come from the same template, sam spec guarantees this
        // read more: https://samtools.github.io/hts-specs/SAMv1.pdf
        r1_names.sort();
        r2_names.sort();
        assert_eq!(r1_names, r2_names, "Mismatch!");
    }

    #[test]
    fn test_paired_end_retention_fetch() {
        let input_path = PathBuf::from("tests/cases/test.paired.bam");

        let target_depth: u32 = 2;
        let output = NamedTempFile::new().unwrap();

        // subsample to 2x depth
        let mut align = Alignment {
            aln: input_path,
            output: Some(output.path().to_path_buf()),
            output_format: Some(Format::Bam),
            coverage: target_depth,
            seed: Some(2109),
            strategy: SubsamplingStrategy::Fetch,
            swap_distance: 5,
            step_size: 100,
            batch_size: 10000,
        };

        align.run().expect("Subsampling failed");

        // verify
        let mut reader = alignment::io::reader::Builder::default()
            .build_from_path(output.path())
            .unwrap();
        // read the header, get the length of chromosome
        let header = reader.read_header().unwrap();

        let mut r1_names: Vec<String> = Vec::new();
        let mut r2_names: Vec<String> = Vec::new();

        for result in reader.records(&header) {
            let record = result.unwrap();
            let name = record.name().unwrap().to_string();
            let flags = record.flags().unwrap();

            if flags.is_first_segment() {
                r1_names.push(name);
            } else if flags.is_last_segment() {
                r2_names.push(name);
            }
        }
        // make sure r1 and r2 names are identical because they come from the same template, sam spec guarantees this
        // read more: https://samtools.github.io/hts-specs/SAMv1.pdf
        r1_names.sort();
        r2_names.sort();

        assert_eq!(r1_names, r2_names, "Mismatch!");
    }
}
