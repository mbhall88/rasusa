mod cli;
mod fastx;

pub use crate::cli::Cli;
pub use crate::fastx::{Fastx, FileType};
use exitfailure::ExitFailure;
use std::io::stdout;
use structopt::StructOpt;

fn main() -> Result<(), ExitFailure> {
    let args = Cli::from_args();

    let input_fastx = Fastx::from_path(&args.input)?;

    let _output_file_handle = match args.output {
        Some(path) => Fastx::from_path(&path)?.create()?,
        None => Box::new(stdout()),
    };

    let _target_total_bases: u64 = args.genome_size * args.coverage;

    let _read_lengths = input_fastx.read_lengths()?;

    Ok(())
}
