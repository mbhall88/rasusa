mod cli;
mod fastx;
mod subsampler;

pub use crate::cli::Cli;
pub use crate::fastx::{Fastx, FileType};
pub use crate::subsampler::SubSampler;

use exitfailure::ExitFailure;
use fern::colors::{Color, ColoredLevelConfig};
use log::{debug, info, warn};
use std::io::stdout;
use structopt::StructOpt;

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

    setup_logger(args.verbose)?;

    debug!("{:?}", args);

    let input_fastx = Fastx::from_path(&args.input[0])?;

    let mut output_file_handle = match args.output.len() {
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

    let target_total_bases: u64 = args.genome_size * args.coverage;
    info!(
        "Target number of bases to subsample to is: {}",
        target_total_bases
    );

    info!("Gathering read lengths...");
    let read_lengths = input_fastx.read_lengths()?;
    info!("{} reads detected", read_lengths.len());

    let subsampler = SubSampler {
        target_total_bases,
        seed: args.seed,
    };

    let mut reads_to_keep = subsampler.indices(&read_lengths);
    info!("Keeping {} reads", reads_to_keep.len());
    debug!("Indices of reads being kept:\n{:?}", reads_to_keep);

    let total_kept_bases: u64 = reads_to_keep
        .iter()
        .map(|&i| read_lengths[i as usize])
        .fold(0, |acc, x| acc + u64::from(x));
    let actual_covg = total_kept_bases / args.genome_size;
    info!("Actual coverage of reads being kept is {:.2}x", actual_covg);

    if let Err(err) = input_fastx.filter_reads_into(&mut reads_to_keep, &mut output_file_handle) {
        warn!("{:?}", err);
    }

    info!("Done.");

    Ok(())
}
