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
    let args = Cli::from_args();

    setup_logger(args.verbose)?;

    debug!("{:?}", args);

    let input_fastx = Fastx::from_path(&args.input)?;

    let mut output_file_handle = match args.output {
        Some(path) => {
            let out_fastx = Fastx::from_path(&path)?;
            if out_fastx.filetype != input_fastx.filetype {
                warn!(
                    "Input ({:?}) and output ({:?}) file types are not the same. \
                     Output will be of the same type as the input.",
                    input_fastx.filetype, out_fastx.filetype
                )
            }
            out_fastx.create()?
        }
        None => Box::new(stdout()),
    };

    let target_total_bases: u64 = args.genome_size * args.coverage;
    info!(
        "Target number of bases to subsample to is: {}",
        target_total_bases
    );

    let read_lengths = input_fastx.read_lengths()?;
    info!("{} reads detected", read_lengths.len());

    let subsampler = SubSampler {
        target_total_bases,
        seed: args.seed,
    };

    let mut reads_to_keep = subsampler.indices(&read_lengths);
    info!("Keeping {} reads", reads_to_keep.len());
    debug!("Indices of reads being kept:\n{:?}", reads_to_keep);

    if let Err(err) = input_fastx.filter_reads_into(&mut reads_to_keep, &mut output_file_handle) {
        warn!("{:?}", err);
    }

    Ok(())
}
