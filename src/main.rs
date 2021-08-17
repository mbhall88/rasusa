use std::io::stdout;

use anyhow::{Context, Result};
use fern::colors::{Color, ColoredLevelConfig};
use log::{debug, error, info};
use structopt::StructOpt;

pub use crate::cli::Cli;
pub use crate::fastx::Fastx;
pub use crate::subsampler::SubSampler;

mod cli;
mod fastx;
mod subsampler;

/// Sets up the logging based on whether you want verbose logging or not. If `verbose` is `false`
/// then info, warning, and error messages will be printed. If `verbose` is `true` then debug
/// messages will also be printed.
///
/// # Errors
/// Will error if `fern` fails to apply the logging setup.
fn setup_logger(verbose: bool) -> Result<(), fern::InitError> {
    let colors = ColoredLevelConfig::new()
        .warn(Color::Yellow)
        .debug(Color::Magenta)
        .error(Color::Red)
        .trace(Color::Green);

    let log_level = if verbose {
        log::LevelFilter::Debug
    } else {
        log::LevelFilter::Info
    };

    fern::Dispatch::new()
        .format(move |out, message, record| {
            out.finish(format_args!(
                "{}[{}][{}] {}",
                chrono::Local::now().format("[%Y-%m-%d][%H:%M:%S]"),
                record.target(),
                colors.color(record.level()),
                message
            ))
        })
        .level(log_level)
        .chain(std::io::stderr())
        .apply()?;
    Ok(())
}

fn main() -> Result<()> {
    let args: Cli = Cli::from_args();
    args.validate_input_output_combination()?;
    let is_illumina = args.input.len() == 2;
    if is_illumina {
        info!("Two input files given. Assuming paired Illumina...")
    }

    setup_logger(args.verbose).context("Failed to setup the logger")?;

    debug!("{:?}", args);

    let input_fastx = Fastx::from_path(&args.input[0]);

    let mut output_handle = match args.output.len() {
        0 => Box::new(stdout()),
        _ => {
            let out_fastx = Fastx::from_path(&args.output[0]);
            out_fastx
                .create()
                .context("unable to create the first output file")?
        }
    };

    let target_total_bases: u64 = args.genome_size * args.coverage;
    info!(
        "Target number of bases to subsample to is: {}",
        target_total_bases
    );

    info!("Gathering read lengths...");
    let mut read_lengths = input_fastx
        .read_lengths()
        .context("unable to gather read lengths for the first input file")?;

    if is_illumina {
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

    let subsampler = SubSampler {
        target_total_bases,
        seed: args.seed,
    };

    let (reads_to_keep, nb_reads_keep) = subsampler.indices(&read_lengths);
    if is_illumina {
        info!("Keeping {} reads from each input", nb_reads_keep);
    } else {
        info!("Keeping {} reads", nb_reads_keep);
    }
    debug!("Indices of reads being kept:\n{:?}", reads_to_keep);

    let mut total_kept_bases =
        input_fastx.filter_reads_into(&reads_to_keep, nb_reads_keep, &mut output_handle)? as u64;

    // repeat the same process for the second input fastx (if illumina)
    if is_illumina {
        let second_input_fastx = Fastx::from_path(&args.input[1]);
        let second_out_fastx = Fastx::from_path(&args.output[1]);
        let mut second_output_handle = second_out_fastx
            .create()
            .context("unable to create the second output file")?;

        total_kept_bases += second_input_fastx.filter_reads_into(
            &reads_to_keep,
            nb_reads_keep,
            &mut second_output_handle,
        )? as u64;
    }

    let actual_covg = total_kept_bases / args.genome_size;
    info!("Actual coverage of kept reads is {:.2}x", actual_covg);

    info!("Done 🎉");

    Ok(())
}
