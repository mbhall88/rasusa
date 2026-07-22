use crate::cli::{
    check_path_exists, parse_compression_format, parse_fraction, parse_level, CliError, Coverage,
    GenomeSize,
};
use crate::fastx::create_output_writer;
use crate::format::{
    default_compression_level, infer_format_from_path, is_fasta_output, output_alignment_format,
    OutputEncoding,
};
use crate::source::{determine_record_source, RecordSource};
use crate::{Runner, SubSampler, SubsampleMode};
use anyhow::{Context, Result};
use clap::Parser;
use log::{debug, info, warn};
use std::io::{stdout, BufWriter};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct Reads {
    /// The FASTA/FASTQ or unaligned SAM/BAM/CRAM file(s) to subsample.
    ///
    /// For paired Illumina, the order matters. i.e., R1 then R2.
    /// Single-file paired-end is also supported for unaligned SAM/BAM/CRAM.
    #[arg(
    value_parser = check_path_exists,
    num_args = 1..=2,
    required = true,
    name = "FILE(S)"
    )]
    pub input: Vec<PathBuf>,

    /// Output filepath(s); stdout if not present.
    ///
    /// For paired Illumina pass this flag twice `-o o1.fq -o o2.fq`  
    ///
    /// NOTE: The order of the pairs is assumed to be the same as the input - e.g., R1 then R2.  
    ///
    /// This option is required for paired input.
    #[arg(short = 'o', long = "output", action = clap::ArgAction::Append)]
    pub output: Vec<PathBuf>,

    /// Genome size to calculate coverage with respect to. e.g., 4.3kb, 7Tb, 9000, 4.1MB
    ///
    /// Alternatively, a FASTA/Q index file can be provided and the genome size will be
    /// set to the sum of all reference sequences.
    ///
    /// If --bases is not provided, this option and --coverage are required
    #[clap(
    short,
    long,
    required_unless_present_any = &["bases", "num", "frac"],
    requires = "coverage",
    value_name = "size|faidx",
    conflicts_with_all = &["num", "frac"]
    )]
    pub genome_size: Option<GenomeSize>,

    /// The desired depth of coverage to subsample the reads to
    ///
    /// If --bases is not provided, this option and --genome-size are required
    #[clap(
    short,
    long,
    value_name = "FLOAT",
    required_unless_present_any = &["bases", "num", "frac"],
    requires = "genome_size",
    conflicts_with_all = &["num", "frac"]
    )]
    pub coverage: Option<Coverage>,

    /// Explicitly set the number of bases required e.g., 4.3kb, 7Tb, 9000, 4.1MB
    ///
    /// If this option is given, --coverage and --genome-size are ignored
    #[clap(short, long, value_name = "bases", conflicts_with_all = &["num", "frac"])]
    pub bases: Option<GenomeSize>,

    /// Subsample to a specific number of reads
    ///
    /// If paired-end reads are passed, this is the number of (matched) reads from EACH file.
    /// This option accepts the same format as genome size - e.g., 1k will take 1000 reads
    #[clap(short, long, value_name = "INT", conflicts_with = "frac")]
    pub num: Option<GenomeSize>,

    /// Subsample to a fraction of the reads - e.g., 0.5 samples half the reads
    ///
    /// Values >1 and <=100 will be automatically converted - e.g., 25 => 0.25
    #[clap(short, long, value_name = "FLOAT", value_parser = parse_fraction, conflicts_with = "num")]
    pub frac: Option<f32>,

    /// Exit with an error if the requested coverage/bases/reads is not possible
    #[clap(short = 'e', long)]
    pub strict: bool,

    /// Random seed to use.
    #[clap(short = 's', long = "seed", value_name = "INT")]
    pub seed: Option<u64>,

    /// Switch on verbosity.
    #[clap(short)]
    pub verbose: bool,

    /// u: uncompressed; b: Bzip2; g: Gzip; l: Lzma; x: Xz (Lzma); z: Zstd
    ///
    /// Rasusa will attempt to infer the output compression format automatically from the filename
    /// extension. This option is used to override that. If writing to stdout, the default is
    /// uncompressed. Note: this is only used for FASTA/FASTQ output.
    #[clap(short = 'Z', long = "compress-type", value_name = "u|b|g|l|x|z", value_parser = parse_compression_format)]
    pub compress_type: Option<niffler::compression::Format>,

    /// Explicitly set the output format.
    ///
    /// If not provided, Rasusa will attempt to infer the format from the filename extension.
    #[clap(short = 'O', long = "output-format", value_enum)]
    pub output_format: Option<crate::cli::OutputFormat>,

    /// Compression level to use if compressing output. Uses the default level for the format if
    /// not specified.
    #[clap(short = 'l', long, value_parser = parse_level, value_name = "1-21")]
    pub compress_level: Option<niffler::Level>,
}

impl Reads {
    /// Checks there is a valid and equal number of `--input` and `--output` arguments given.
    ///
    /// # Errors
    /// A [`CliError::BadInputOutputCombination`](#clierror) is returned for the following:
    /// - Either `--input` or `--output` are passed more than twice
    /// - An unequal number of `--input` and `--output` are passed. The only exception to
    ///   this is if one `--input` and zero `--output` are passed, in which case, the output
    ///   will be sent to STDOUT.
    pub fn validate_input_output_combination(&self) -> std::result::Result<(), CliError> {
        let out_len = self.output.len();
        let in_len = self.input.len();

        if in_len > 2 {
            let msg = String::from("Got more than 2 files for input.");
            return Err(CliError::BadInputOutputCombination(msg));
        }
        if out_len > 2 {
            let msg = String::from("Got more than 2 files for output.");
            return Err(CliError::BadInputOutputCombination(msg));
        }
        match in_len as isize - out_len as isize {
            diff if diff == 1 && in_len == 1 => Ok(()),
            diff if diff != 0 => Err(CliError::BadInputOutputCombination(format!(
                "Got {in_len} --input but {out_len} --output"
            ))),
            _ => Ok(()),
        }
    }
}

impl Runner for Reads {
    fn run(&mut self) -> Result<()> {
        self.validate_input_output_combination()?;
        let is_paired = self.input.len() == 2;
        if is_paired {
            info!("Two input files given. Assuming paired Illumina...")
        }

        let input_source = determine_record_source(&self.input[0]);
        let input_format = infer_format_from_path(&self.input[0]);

        let check_conversion = |in_fmt: Option<noodles_util::alignment::io::Format>,
                                out_path: Option<&std::path::PathBuf>|
         -> Result<()> {
            if in_fmt.is_none() {
                let explicit_out_fmt = output_alignment_format(self.output_format.as_ref());
                let inferred_out_fmt = out_path.and_then(|p| infer_format_from_path(p));
                let has_alignment_output = explicit_out_fmt.is_some() || inferred_out_fmt.is_some();

                if has_alignment_output {
                    let out_fmt = explicit_out_fmt.or(inferred_out_fmt).unwrap();
                    return Err(anyhow::anyhow!(
                        "Conversion from FASTA/FASTQ to {:?} is not supported. Please use FASTA/FASTQ output for FASTA/FASTQ input.",
                        out_fmt
                    ));
                }
            }
            Ok(())
        };

        check_conversion(input_format, self.output.first())?;
        if is_paired {
            let second_input_format = infer_format_from_path(&self.input[1]);
            check_conversion(second_input_format, self.output.get(1))?;
        }

        let mut output_handle = match self.output.len() {
            0 => match self.compress_type {
                None => Box::new(BufWriter::new(stdout().lock())) as Box<dyn std::io::Write>,
                Some(fmt) => {
                    let lvl = default_compression_level(fmt);
                    niffler::basic::get_writer(Box::new(BufWriter::new(stdout().lock())), fmt, lvl)?
                }
            },
            _ => create_output_writer(&self.output[0], self.compress_level, self.compress_type)
                .context("unable to create the first output file")?,
        };

        let second_input_source = if is_paired {
            Some(determine_record_source(&self.input[1]))
        } else {
            None
        };

        self.subsample(
            input_source.as_ref(),
            second_input_source.as_deref(),
            &mut *output_handle,
        )
    }
}

impl Reads {
    /// Runs the core subsampling orchestration against already-constructed sources and sinks.
    ///
    /// This is the in-process seam used by both [`Runner::run`] (which builds the first
    /// source/sink from `self.input[0]`/`self.output`) and unit tests (which can inject a source
    /// backed by a temp file and an in-memory `Vec<u8>` sink, without spawning the CLI as a
    /// subprocess). The second (paired-mode) output is always file-backed - as it was before this
    /// seam existed - so it's still constructed from `self.output[1]` here rather than injected.
    fn subsample(
        &self,
        input_source: &dyn RecordSource,
        second_input_source: Option<&dyn RecordSource>,
        output_handle: &mut dyn std::io::Write,
    ) -> Result<()> {
        let is_paired = second_input_source.is_some();
        let input_format = infer_format_from_path(&self.input[0]);

        let target_total_bases: Option<u64> = match (self.genome_size, self.coverage, self.bases) {
            (_, _, Some(bases)) => Some(u64::from(bases)),
            (Some(gsize), Some(cov), _) => Some(gsize * cov),
            _ => None,
        };

        if let Some(bases) = target_total_bases {
            info!("Target number of bases to subsample to is: {bases}",);
        }

        let (mode, total_reads, read_lengths) = if let Some(ttb) = target_total_bases {
            info!("Gathering read lengths...");
            let mut lengths = input_source
                .read_lengths()
                .context("unable to gather read lengths for the first input file")?;

            if let Some(second_source) = second_input_source {
                let expected_num_reads = lengths.len();
                info!("Gathering read lengths for second input file...");
                let mate_lengths = second_source
                    .read_lengths()
                    .context("unable to gather read lengths for the second input file")?;

                check_paired_counts(expected_num_reads, mate_lengths.len())?;
                // add the paired read lengths to the existing lengths
                for (i, len) in mate_lengths.iter().enumerate() {
                    lengths[i] += len;
                }
            }
            info!("{} reads detected", lengths.len());

            let total_input_bases: u64 = lengths.iter().map(|&x| x as u64).sum();

            // calculate the depth of coverage if using coverage-based subsampling
            if let Some(size) = self.genome_size {
                let depth_of_covg = (total_input_bases as f64) / f64::from(size);
                info!("Input coverage is {:.2}x", depth_of_covg);

                if self.strict && Coverage(depth_of_covg as f32) < self.coverage.unwrap() {
                    return Err(anyhow::anyhow!(
                        "Requested coverage ({:.2}x) is not possible as the actual coverage is {:.2}x",
                        self.coverage.unwrap().0,
                        depth_of_covg
                    ));
                }
            }

            if let Some(req_bases) = self.bases {
                let req_bases = u64::from(req_bases);
                if self.strict && req_bases > total_input_bases {
                    return Err(anyhow::anyhow!(
                        "Requested number of bases ({}) is more than the input ({})",
                        req_bases,
                        total_input_bases
                    ));
                }
            }

            let total_reads = lengths.len();
            (SubsampleMode::ByBases(ttb), total_reads, lengths)
        } else {
            // in this mode, selection only depends on the read *count* - avoid gathering the
            // (unused) per-read lengths vector.
            info!("Counting reads...");
            let first_count = input_source
                .count()
                .context("unable to count reads in the first input file")?;

            let total_reads = if let Some(second_source) = second_input_source {
                info!("Counting reads in second input file...");
                let second_count = second_source
                    .count()
                    .context("unable to count reads in the second input file")?;

                check_paired_counts(first_count, second_count)?;
                first_count
            } else {
                first_count
            };
            info!("{} reads detected", total_reads);

            let n = match (self.num, self.frac) {
                (Some(n), None) => u64::from(n),
                (None, Some(f)) => {
                    let n = ((f as f64) * (total_reads as f64)).round() as u64;
                    if n == 0 {
                        if !self.strict {
                            warn!(
                                "Requested fraction of reads ({} * {}) was rounded to 0",
                                f, total_reads
                            );
                        } else {
                            return Err(anyhow::anyhow!(
                                "Requested fraction of reads ({} * {}) was rounded to 0",
                                f,
                                total_reads
                            ));
                        }
                    }
                    n
                }
                _ => {
                    return Err(anyhow::anyhow!(
                        "Either --num or --frac must be given when not subsampling by --bases/--coverage"
                    ));
                }
            };

            if self.strict && n as usize > total_reads {
                return Err(anyhow::anyhow!(
                    "Requested number of reads ({}) is more than the input ({})",
                    n,
                    total_reads
                ));
            }
            info!("Target number of reads to subsample to is: {}", n);

            (SubsampleMode::ByReads(n), total_reads, Vec::new())
        };

        let subsampler = SubSampler {
            mode,
            seed: self.seed,
        };

        let (reads_to_keep, nb_reads_to_keep) = subsampler.indices(total_reads, &read_lengths);
        if is_paired {
            info!("Keeping {} reads from each input", nb_reads_to_keep);
        } else {
            info!("Keeping {} reads", nb_reads_to_keep);
        }
        const DEBUG_MAX_LOGGED_INDICES: usize = 10_000;
        if reads_to_keep.len() <= DEBUG_MAX_LOGGED_INDICES {
            debug!("Indices of reads being kept:\n{:?}", reads_to_keep);
        } else {
            debug!(
                "Indices of reads being kept: omitted ({} reads exceeds debug print cap of {})",
                reads_to_keep.len(),
                DEBUG_MAX_LOGGED_INDICES
            );
        }

        let output_format_1 = output_alignment_format(self.output_format.as_ref()).or_else(|| {
            if self.output.is_empty() {
                input_format
            } else {
                infer_format_from_path(&self.output[0])
            }
        });

        let encoding_1 = match output_format_1 {
            Some(fmt) => OutputEncoding::Alignment(fmt),
            None => {
                let fallback_path = if self.output.is_empty() {
                    &self.input[0]
                } else {
                    &self.output[0]
                };
                OutputEncoding::Fastx {
                    fasta: is_fasta_output(self.output_format.as_ref(), fallback_path),
                }
            }
        };

        let mut total_kept_bases = input_source.filter_reads_into(
            &reads_to_keep,
            nb_reads_to_keep,
            output_handle,
            encoding_1,
        )? as u64;

        // repeat the same process for the second input (if illumina)
        if let Some(second_input_source) = second_input_source {
            let second_input_format = infer_format_from_path(&self.input[1]);

            let mut second_output_handle =
                create_output_writer(&self.output[1], self.compress_level, self.compress_type)
                    .context("unable to create the second output file")?;

            let output_format_2 =
                output_alignment_format(self.output_format.as_ref()).or_else(|| {
                    if self.output.len() < 2 {
                        second_input_format
                    } else {
                        infer_format_from_path(&self.output[1])
                    }
                });

            let encoding_2 = match output_format_2 {
                Some(fmt) => OutputEncoding::Alignment(fmt),
                None => {
                    let fallback_path = if self.output.len() < 2 {
                        &self.input[1]
                    } else {
                        &self.output[1]
                    };
                    OutputEncoding::Fastx {
                        fasta: is_fasta_output(self.output_format.as_ref(), fallback_path),
                    }
                }
            };

            total_kept_bases += second_input_source.filter_reads_into(
                &reads_to_keep,
                nb_reads_to_keep,
                &mut second_output_handle,
                encoding_2,
            )? as u64;
        }

        if let Some(gsize) = self.genome_size {
            let actual_covg = total_kept_bases / gsize;
            // safe to unwrap self.coverage as we have already ensured genome size and coverage co-occur
            if Coverage(actual_covg as f32) < self.coverage.unwrap() {
                warn!(
                "Requested coverage ({:.2}x) is not possible as the actual coverage is {:.2}x - \
                output will be the same as the input",
                self.coverage.unwrap().0,
                actual_covg
            );
            } else {
                info!("Actual coverage of kept reads is {:.2}x", actual_covg);
            }
        } else {
            info!("Kept {} bases", total_kept_bases);
        }

        info!("Done 🎉");
        Ok(())
    }
}

/// Checks that a paired input's two files report the same number of reads, logging a success
/// message on match. Both "by bases" (per-file `read_lengths().len()`) and "by reads"
/// (per-file `count()`) branches of [`Reads::subsample`] hit the same mismatch condition and
/// want the same error/log wording.
fn check_paired_counts(first_count: usize, second_count: usize) -> Result<()> {
    if second_count != first_count {
        return Err(anyhow::anyhow!(
            "First input has {} reads, but the second has {} reads. Paired Illumina files are assumed to have the same number of reads.",
            first_count,
            second_count
        ));
    }
    info!(
        "Both input files have the same number of reads ({}) 👍",
        first_count
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fastx::Fastx;
    use assert_cmd::Command;
    use std::io::Write as _;
    use std::str::FromStr;
    use tempfile::NamedTempFile;

    const SUB: &str = "reads";

    /// Builds a temp FASTQ file with `n` reads, each `len` bases long, for use with the
    /// in-process `subsample` seam (mirrors the pattern `alignment.rs` uses for its own
    /// in-process tests, but for the reads path).
    fn write_fastq(n: usize, len: usize) -> NamedTempFile {
        let mut file = NamedTempFile::with_suffix(".fq").unwrap();
        let seq = "A".repeat(len);
        let qual = "I".repeat(len);
        for i in 0..n {
            writeln!(file, "@read{i}\n{seq}\n+\n{qual}").unwrap();
        }
        file.flush().unwrap();
        file
    }

    fn default_reads(input: Vec<PathBuf>) -> Reads {
        Reads {
            input,
            output: vec![],
            genome_size: None,
            coverage: None,
            bases: None,
            num: None,
            frac: None,
            strict: false,
            seed: Some(1),
            verbose: false,
            compress_type: None,
            output_format: None,
            compress_level: None,
        }
    }

    #[test]
    fn frac_rounds_to_zero_without_strict_warns_and_keeps_nothing() {
        let input = write_fastq(2, 10);
        let source = Fastx::from_path(input.path());
        let mut reads = default_reads(vec![input.path().to_path_buf()]);
        reads.frac = Some(0.01); // 0.01 * 2 rounds to 0

        let mut output = Vec::new();
        let result = reads.subsample(&source, None, &mut output);

        assert!(result.is_ok());
        assert!(output.is_empty());
    }

    #[test]
    fn frac_rounds_to_zero_with_strict_returns_err() {
        let input = write_fastq(2, 10);
        let source = Fastx::from_path(input.path());
        let mut reads = default_reads(vec![input.path().to_path_buf()]);
        reads.frac = Some(0.01);
        reads.strict = true;

        let mut output = Vec::new();
        let result = reads.subsample(&source, None, &mut output);

        assert!(result.is_err());
    }

    #[test]
    fn num_more_than_available_with_strict_returns_err() {
        let input = write_fastq(2, 10);
        let source = Fastx::from_path(input.path());
        let mut reads = default_reads(vec![input.path().to_path_buf()]);
        reads.num = Some(GenomeSize::from_str("10").unwrap());
        reads.strict = true;

        let mut output = Vec::new();
        let result = reads.subsample(&source, None, &mut output);

        assert!(result.is_err());
    }

    #[test]
    fn num_more_than_available_without_strict_keeps_all() {
        let input = write_fastq(2, 10);
        let source = Fastx::from_path(input.path());
        let mut reads = default_reads(vec![input.path().to_path_buf()]);
        reads.num = Some(GenomeSize::from_str("10").unwrap());

        let mut output = Vec::new();
        let result = reads.subsample(&source, None, &mut output);

        assert!(result.is_ok());
        assert_eq!(output.iter().filter(|&&b| b == b'@').count(), 2);
    }

    #[test]
    fn insufficient_coverage_without_strict_warns_instead_of_erroring() {
        // 2 reads * 10 bases = 20 bases total, requesting 5000x coverage of a 1000bp genome
        // (5,000,000 bases) is unattainable - this should warn, not error, when not strict.
        let input = write_fastq(2, 10);
        let source = Fastx::from_path(input.path());
        let mut reads = default_reads(vec![input.path().to_path_buf()]);
        reads.genome_size = Some(GenomeSize::from_str("1000").unwrap());
        reads.coverage = Some(Coverage(5000.0));

        let mut output = Vec::new();
        let result = reads.subsample(&source, None, &mut output);

        // all reads get kept (can never reach the target), but this is not an error
        assert!(result.is_ok());
        assert_eq!(output.iter().filter(|&&b| b == b'@').count(), 2);
    }

    #[test]
    fn insufficient_coverage_with_strict_returns_err() {
        let input = write_fastq(2, 10);
        let source = Fastx::from_path(input.path());
        let mut reads = default_reads(vec![input.path().to_path_buf()]);
        reads.genome_size = Some(GenomeSize::from_str("1000").unwrap());
        reads.coverage = Some(Coverage(5000.0));
        reads.strict = true;

        let mut output = Vec::new();
        let result = reads.subsample(&source, None, &mut output);

        assert!(result.is_err());
    }

    #[test]
    fn paired_read_count_mismatch_returns_err_instead_of_exiting() {
        let input1 = write_fastq(3, 10);
        let input2 = write_fastq(2, 10);
        let source1 = Fastx::from_path(input1.path());
        let source2 = Fastx::from_path(input2.path());
        let mut reads = default_reads(vec![
            input1.path().to_path_buf(),
            input2.path().to_path_buf(),
        ]);
        reads.num = Some(GenomeSize::from_str("1").unwrap());

        let mut output = Vec::new();
        let result = reads.subsample(&source1, Some(&source2), &mut output);

        let err = result.expect_err("mismatched paired read counts should be an error");
        assert!(err.to_string().contains("First input has"));
    }

    #[test]
    fn paired_read_count_mismatch_by_bases_returns_err_instead_of_exiting() {
        let input1 = write_fastq(3, 10);
        let input2 = write_fastq(2, 10);
        let source1 = Fastx::from_path(input1.path());
        let source2 = Fastx::from_path(input2.path());
        let mut reads = default_reads(vec![
            input1.path().to_path_buf(),
            input2.path().to_path_buf(),
        ]);
        reads.bases = Some(GenomeSize::from_str("10").unwrap());

        let mut output = Vec::new();
        let result = reads.subsample(&source1, Some(&source2), &mut output);

        let err = result.expect_err("mismatched paired read counts should be an error");
        assert!(err.to_string().contains("First input has"));
    }

    #[test]
    fn no_args_given_raises_error() {
        let passed_args = vec![SUB];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn no_input_file_given_raises_error() {
        let passed_args = vec![SUB, "-c", "30", "-g", "3mb"];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn no_coverage_given_raises_error() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![SUB, infile, "-g", "3mb"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn invalid_coverage_given_raises_error() {
        let passed_args = vec![SUB, "in.fq", "-g", "3mb", "-c", "foo"];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn no_genome_size_given_raises_error() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![SUB, infile, "-c", "5"];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn no_genome_size_but_bases_and_coverage() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![SUB, infile, "-b", "5", "-c", "7"];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn bases_and_coverage_and_genome_size_all_allowed() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![SUB, infile, "-b", "5", "-c", "7", "-g", "5m"];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().success();
    }

    #[test]
    fn no_genome_size_or_coverage_given_but_bases_given() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![SUB, infile, "-b", "5"];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().success();
    }

    #[test]
    fn no_genome_size_or_coverage_given_but_num_given() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![SUB, infile, "-n", "5"];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().success();
    }

    #[test]
    fn no_genome_size_or_coverage_given_but_frac_given() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![SUB, infile, "-f", "5"];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().success();
    }

    #[test]
    fn num_and_coverage_and_genome_size_not_allowed() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![SUB, infile, "-n", "5", "-c", "7", "-g", "5m"];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn frac_and_coverage_and_genome_size_not_allowed() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![SUB, infile, "-f", "5", "-c", "7", "-g", "5m"];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn frac_and_bases_not_allowed() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![SUB, infile, "-f", "5", "-b", "7"];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn frac_and_num_not_allowed() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![SUB, infile, "-f", "5", "-n", "7"];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn bases_and_num_not_allowed() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![SUB, infile, "-b", "5", "-n", "7"];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn invalid_genome_size_given_raises_error() {
        let passed_args = vec![SUB, "in.fq", "-c", "5", "-g", "8jb"];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn faidx_given_as_genome_size() {
        let infile = "tests/cases/r1.fq.gz";
        let faidx = "tests/cases/h37rv.fa.fai";
        let passed_args = vec![SUB, infile, "-c", "5", "-g", faidx];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().success();
    }

    #[test]
    fn invalid_seed_given_raises_error() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![SUB, infile, "-c", "5", "-g", "8mb", "-s", "foo"];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn all_valid_args_parsed_as_expected() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![
            SUB,
            infile,
            "-c",
            "5",
            "-g",
            "8mb",
            "-s",
            "88",
            "-o",
            "/dev/null",
        ];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().success();
    }

    #[test]
    fn two_outputs_passed_with_one_option_flag() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![
            SUB,
            infile,
            infile,
            "-c",
            "5",
            "-g",
            "8mb",
            "-s",
            "88",
            "-o",
            "/dev/null",
            "/dev/null",
        ];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn three_inputs_raises_error() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![
            SUB,
            infile,
            infile,
            infile,
            "-c",
            "5",
            "-g",
            "8mb",
            "-s",
            "88",
            "-o",
            "my/output/file.fq",
        ];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn three_outputs_raises_error() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![
            SUB,
            infile,
            infile,
            "-c",
            "5",
            "-g",
            "8mb",
            "-s",
            "88",
            "-o",
            "my/output/file.fq",
            "-o",
            "out.fq",
            "-o",
            "out.fq",
        ];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn one_input_no_output_is_ok() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![SUB, infile, "-c", "5", "-g", "8mb", "-s", "88"];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().success();
    }

    #[test]
    fn two_inputs_one_output_raises_error() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![
            SUB, infile, infile, "-c", "5", "-g", "8mb", "-s", "88", "-o", "out.fq",
        ];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn one_input_two_outputs_raises_error() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![
            SUB, infile, "-c", "5", "-g", "8mb", "-s", "88", "-o", "out.fq", "-o", "out.fq",
        ];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn two_input_two_outputs_is_ok() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![
            SUB, infile, infile, "-c", "5", "-g", "8mb", "-s", "88", "-o", "out.fq", "-o", "out.fq",
        ];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().success();
    }

    #[test]
    fn two_input_two_outputs_is_ok_when_positional_args_at_end() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![
            SUB, "-c", "5", "-g", "8mb", "-s", "88", "-o", "out.fq", "-o", "out.fq", infile, infile,
        ];

        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
        cmd.args(passed_args).assert().success();
    }
}
