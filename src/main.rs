mod cli;
mod fastx;

pub use crate::cli::Cli;
pub use crate::fastx::{Fastx, FileType};
use exitfailure::ExitFailure;
use std::io::stdout;
use structopt::StructOpt;

fn main() -> Result<(), ExitFailure> {
    let args = Cli::from_args();
    // determine if file is fasta or fastq
    let input_fastx = Fastx::from_path(&args.input)?;
    // get file handle for output file/stdout
    let _output_file_handle = match args.output {
        Some(path) => Fastx::from_path(&path)?.create()?,
        None => Box::new(stdout()),
    };
    // get input file read lengths
    let _read_lengths = input_fastx.read_lengths()?;

    Ok(())
}
