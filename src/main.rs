#![allow(clippy::redundant_clone)]

extern crate core;

// due to a structopt problem
use std::io::stdout;

use anyhow::{Context, Result};
use clap::Parser;
use env_logger::Builder;
use log::{debug, error, info, warn, LevelFilter};
use niffler::compression;

pub use crate::cli::Cli;
use crate::cli::Coverage;
pub use crate::fastx::Fastx;
pub use crate::subsampler::SubSampler;

mod cli;
mod fastx;
mod subsampler;

fn main() -> Result<()> {
    let args: Cli = Cli::parse();
    // Initialize logger
    let log_lvl = if args.verbose {
        LevelFilter::Debug
    } else {
        LevelFilter::Info
    };
    let mut log_builder = Builder::new();
    log_builder
        .filter(None, log_lvl)
        .format_module_path(false)
        .format_target(false)
        .init();

    debug!("{:?}", args);

    args.validate_input_output_combination()?;
    let is_paired = args.input.len() == 2;
    if is_paired {
        info!("Two input files given. Assuming paired Illumina...")
    }

    let input_fastx = Fastx::from_path(&args.input[0]);

    let mut output_handle = match args.output.len() {
        0 => match args.output_type {
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
            let out_fastx = Fastx::from_path(&args.output[0]);
            out_fastx
                .create(args.compress_level, args.output_type)
                .context("unable to create the first output file")?
        }
    };

    let target_total_bases: Option<u64> = match (args.genome_size, args.coverage, args.bases) {
        (_, _, Some(bases)) => Some(u64::from(bases)),
        (Some(gsize), Some(cov), _) => Some(gsize * cov),
        _ => None,
    };

    if target_total_bases.is_some() {
        info!(
            "Target number of bases to subsample to is: {}",
            target_total_bases.unwrap()
        );
    }

    info!("Gathering read lengths...");
    let mut read_lengths = input_fastx
        .read_lengths()
        .context("unable to gather read lengths for the first input file")?;

    if is_paired {
        let second_input_fastx = Fastx::from_path(&args.input[1]);
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
    if args.genome_size.is_some() {
        let number_of_bases: u64 = read_lengths.iter().map(|&x| x as u64).sum();
        let depth_of_covg = (number_of_bases as f64) / f64::from(args.genome_size.unwrap());
        info!("Input coverage is {:.2}x", depth_of_covg);
    }

    let num_reads = match (args.num, args.frac) {
        (Some(n), None) => Some(u64::from(n)),
        (None, Some(f)) => {
            let n = ((f as f64) * (read_lengths.len() as f64)).round() as u64;
            if n == 0 {
                warn!(
                    "Requested fraction of reads ({} * {}) was rounded to 0",
                    f,
                    read_lengths.len()
                );
            }
            Some(n)
        }
        _ => None,
    };

    let subsampler = SubSampler {
        target_total_bases,
        seed: args.seed,
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
        input_fastx.filter_reads_into(&reads_to_keep, nb_reads_to_keep, &mut output_handle)? as u64;

    // repeat the same process for the second input fastx (if illumina)
    if is_paired {
        let second_input_fastx = Fastx::from_path(&args.input[1]);
        let second_out_fastx = Fastx::from_path(&args.output[1]);
        let mut second_output_handle = second_out_fastx
            .create(args.compress_level, args.output_type)
            .context("unable to create the second output file")?;

        total_kept_bases += second_input_fastx.filter_reads_into(
            &reads_to_keep,
            nb_reads_to_keep,
            &mut second_output_handle,
        )? as u64;
    }

    if let Some(gsize) = args.genome_size {
        let actual_covg = total_kept_bases / gsize;
        // safe to unwrap args.coverage as we have already ensured genome size and coverage co-occur
        if Coverage(actual_covg as f32) < args.coverage.unwrap() {
            warn!(
                "Requested coverage ({:.2}x) is not possible as the actual coverage is {:.2}x - \
                output will be the same as the input",
                args.coverage.unwrap().0,
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
