mod cli;
mod fastx;
mod subsampler;

pub use crate::cli::Cli;
pub use crate::fastx::{Fastx, FileType};
pub use crate::subsampler::SubSampler;

use exitfailure::ExitFailure;
use std::io::stdout;
use structopt::StructOpt;

fn main() -> Result<(), ExitFailure> {
    let args = Cli::from_args();

    let input_fastx = Fastx::from_path(&args.input)?;

    let mut output_file_handle = match args.output {
        Some(path) => Fastx::from_path(&path)?.create()?,
        None => Box::new(stdout()),
    };

    let target_total_bases: u64 = args.genome_size * args.coverage;

    let read_lengths = input_fastx.read_lengths()?;

    let subsampler = SubSampler {
        target_total_bases,
        seed: args.seed,
    };

    let mut reads_to_keep = subsampler.indices(&read_lengths);

    if let Err(err) = input_fastx.filter_reads_into(&mut reads_to_keep, &mut output_file_handle) {
        // todo: add logging warning
        println!("{:?}", err);
    }

    Ok(())
}
