use crate::cli::{
    check_path_exists, parse_compression_format, parse_fraction, parse_level, CliError, Coverage,
    GenomeSize,
};
use crate::{Fastx, Runner, SubSampler};
use anyhow::{Context, Result};
use clap::Parser;
use log::{debug, error, info, warn};
use niffler::compression;
use std::io::stdout;
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct Reads {
    /// The fast{a,q} file(s) to subsample.
    ///
    /// For paired Illumina, the order matters. i.e., R1 then R2.
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
    /// uncompressed
    #[clap(short = 'O', long, value_name = "u|b|g|l|x|z", value_parser = parse_compression_format)]
    pub output_type: Option<niffler::compression::Format>,

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

        let input_fastx = Fastx::from_path(&self.input[0]);

        let mut output_handle = match self.output.len() {
            0 => match self.output_type {
                None => Box::new(stdout()),
                Some(fmt) => {
                    let lvl = match fmt {
                        compression::Format::Gzip => compression::Level::Six,
                        compression::Format::Bzip => compression::Level::Nine,
                        compression::Format::Lzma => compression::Level::Six,
                        compression::Format::Zstd => compression::Level::Three,
                        _ => compression::Level::Zero,
                    };
                    niffler::basic::get_writer(Box::new(stdout()), fmt, lvl)?
                }
            },
            _ => {
                let out_fastx = Fastx::from_path(&self.output[0]);
                out_fastx
                    .create(self.compress_level, self.output_type)
                    .context("unable to create the first output file")?
            }
        };

        let target_total_bases: Option<u64> = match (self.genome_size, self.coverage, self.bases) {
            (_, _, Some(bases)) => Some(u64::from(bases)),
            (Some(gsize), Some(cov), _) => Some(gsize * cov),
            _ => None,
        };

        if let Some(bases) = target_total_bases {
            info!("Target number of bases to subsample to is: {bases}",);
        }

        info!("Gathering read lengths...");
        let mut read_lengths = input_fastx
            .read_lengths()
            .context("unable to gather read lengths for the first input file")?;

        if is_paired {
            let second_input_fastx = Fastx::from_path(&self.input[1]);
            let expected_num_reads = read_lengths.len();
            info!("Gathering read lengths for second input file...");
            let mate_lengths = second_input_fastx
                .read_lengths()
                .context("unable to gather read lengths for the second input file")?;

            if mate_lengths.len() != expected_num_reads {
                error!("First input has {} reads, but the second has {} reads. Paired Illumina files are assumed to have the same number of reads. The results of this subsample may not be as expected now.", expected_num_reads, read_lengths.len());
                std::process::exit(1);
            } else {
                info!(
                    "Both input files have the same number of reads ({}) 👍",
                    expected_num_reads
                );
            }
            // add the paired read lengths to the existing lengths
            for (i, len) in mate_lengths.iter().enumerate() {
                read_lengths[i] += len;
            }
        }
        info!("{} reads detected", read_lengths.len());

        // calculate the depth of coverage if using coverage-based subsampling
        if let Some(size) = self.genome_size {
            let number_of_bases: u64 = read_lengths.iter().map(|&x| x as u64).sum();
            let depth_of_covg = (number_of_bases as f64) / f64::from(size);
            info!("Input coverage is {:.2}x", depth_of_covg);

            if self.strict
                && target_total_bases.is_some()
                && Coverage(depth_of_covg as f32) < self.coverage.unwrap()
            {
                return Err(anyhow::anyhow!(
                    "Requested coverage ({:.2}x) is not possible as the actual coverage is {:.2}x",
                    self.coverage.unwrap().0,
                    depth_of_covg
                ));
            }
        }

        if let Some(req_bases) = self.bases {
            let total_bases: u64 = read_lengths.iter().map(|&x| x as u64).sum();
            let req_bases = u64::from(req_bases);
            if self.strict && req_bases > total_bases {
                return Err(anyhow::anyhow!(
                    "Requested number of bases ({}) is more than the input ({})",
                    req_bases,
                    total_bases
                ));
            }
        }

        let num_reads = match (self.num, self.frac) {
            (Some(n), None) => Some(u64::from(n)),
            (None, Some(f)) => {
                let n = ((f as f64) * (read_lengths.len() as f64)).round() as u64;
                if n == 0 {
                    if !self.strict {
                        warn!(
                            "Requested fraction of reads ({} * {}) was rounded to 0",
                            f,
                            read_lengths.len()
                        );
                    } else {
                        return Err(anyhow::anyhow!(
                            "Requested fraction of reads ({} * {}) was rounded to 0",
                            f,
                            read_lengths.len()
                        ));
                    }
                }
                Some(n)
            }
            _ => None,
        };

        if let Some(n) = num_reads {
            if self.strict && n as usize > read_lengths.len() {
                return Err(anyhow::anyhow!(
                    "Requested number of reads ({}) is more than the input ({})",
                    n,
                    read_lengths.len()
                ));
            }
            info!("Target number of reads to subsample to is: {}", n);
        }

        let subsampler = SubSampler {
            target_total_bases,
            seed: self.seed,
            num_reads,
        };

        let (reads_to_keep, nb_reads_to_keep) = subsampler.indices(&read_lengths);
        if is_paired {
            info!("Keeping {} reads from each input", nb_reads_to_keep);
        } else {
            info!("Keeping {} reads", nb_reads_to_keep);
        }
        debug!("Indices of reads being kept:\n{:?}", reads_to_keep);

        let mut total_kept_bases =
            input_fastx.filter_reads_into(&reads_to_keep, nb_reads_to_keep, &mut output_handle)?
                as u64;

        // repeat the same process for the second input fastx (if illumina)
        if is_paired {
            let second_input_fastx = Fastx::from_path(&self.input[1]);
            let second_out_fastx = Fastx::from_path(&self.output[1]);
            let mut second_output_handle = second_out_fastx
                .create(self.compress_level, self.output_type)
                .context("unable to create the second output file")?;

            total_kept_bases += second_input_fastx.filter_reads_into(
                &reads_to_keep,
                nb_reads_to_keep,
                &mut second_output_handle,
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

#[cfg(test)]
mod tests {
    use assert_cmd::Command;

    const SUB: &str = "reads";

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
