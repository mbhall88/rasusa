use std::io::stdout;

use exitfailure::ExitFailure;
use fern::colors::{Color, ColoredLevelConfig};
use log::{debug, error, info, warn};
use structopt::StructOpt;

pub use crate::cli::Cli;
pub use crate::fastx::{Fastx, FileType};
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

fn main() -> Result<(), ExitFailure> {
    let args: Cli = Cli::from_args();
    args.validate_input_output_combination()?;
    let is_illumina = args.input.len() == 2;
    if is_illumina {
        info!("Two input files given. Assuming paired Illumina...")
    }

    setup_logger(args.verbose)?;

    debug!("{:?}", args);

    let input_fastx = Fastx::from_path(&args.input[0])?;

    let mut output_handle = match args.output.len() {
        0 => Box::new(stdout()),
        _ => {
            let out_fastx = Fastx::from_path(&args.output[0])?;
            if out_fastx.filetype != input_fastx.filetype {
                warn!(
                    "Input ({:?}) and output ({:?}) file types are not the same. \
                     Output will be of the same type as the input.",
                    input_fastx.filetype, out_fastx.filetype
                )
            }
            out_fastx.create()?
        }
    };

    let mut target_total_bases: u64 = args.genome_size * args.coverage;
    info!(
        "Target number of bases to subsample to is: {}",
        target_total_bases
    );

    info!("Gathering read lengths...");
    let mut read_lengths = input_fastx.read_lengths()?;

    if is_illumina {
        info!("{} reads detected in the first input", read_lengths.len());
        target_total_bases /= 2;
    } else {
        info!("{} reads detected", read_lengths.len());
    }

    let subsampler = SubSampler {
        target_total_bases,
        seed: args.seed,
    };

    let reads_to_keep = subsampler.indices(&read_lengths);
    if is_illumina {
        info!("Keeping {} reads from each input", reads_to_keep.len());
    } else {
        info!("Keeping {} reads", reads_to_keep.len());
    }
    debug!("Indices of reads being kept:\n{:?}", reads_to_keep);

    if let Err(err) = input_fastx.filter_reads_into(reads_to_keep.clone(), &mut output_handle) {
        warn!("{:?}", err);
    }

    let mut total_kept_bases: u64 = reads_to_keep
        .iter()
        .map(|&i| read_lengths[i as usize])
        .fold(0, |acc, x| acc + u64::from(x));

    // repeat the same process for the second input fastx (if illumina)
    if is_illumina {
        let second_input_fastx = Fastx::from_path(&args.input[1])?;
        let second_out_fastx = Fastx::from_path(&args.output[1])?;
        let mut second_output_handle = second_out_fastx.create()?;

        let expected_num_reads = read_lengths.len();
        info!("Gathering read lengths for second input file...");
        read_lengths = second_input_fastx.read_lengths()?;

        if read_lengths.len() != expected_num_reads {
            error!("First input has {} reads, but the second has {} reads. Paired Illumina files are assumed to have the same number of reads. The results of this subsample may not be as expected now.", expected_num_reads, read_lengths.len())
        } else {
            info!(
                "Both input files have the same number of reads ({}) 👍",
                expected_num_reads
            )
        }

        total_kept_bases += reads_to_keep
            .iter()
            .map(|&i| read_lengths[i as usize])
            .fold(0, |acc, x| acc + u64::from(x));

        if let Err(err) =
            second_input_fastx.filter_reads_into(reads_to_keep, &mut second_output_handle)
        {
            warn!("{:?}", err);
        }
    }

    let actual_covg = total_kept_bases / args.genome_size;
    info!("Actual coverage of kept reads is {:.2}x", actual_covg);

    info!("Done 🎉");

    Ok(())
}
