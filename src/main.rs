mod cli;
mod fastx;

pub use crate::cli::Cli;
pub use crate::fastx::Fastx;
use exitfailure::ExitFailure;
use std::io::stdout;
use structopt::StructOpt;

fn main() -> Result<(), ExitFailure> {
    let args = Cli::from_args();
    // determine if file is fasta or fastq
    let input_fastx = Fastx::from_path(&args.input)?;
    let _input_file_handle = input_fastx.open()?;
    // get file handle for output file/stdout
    let _output_file_handle = match args.output {
        Some(path) => Fastx::from_path(&path)?.create()?,
        None => Box::new(stdout()),
    };
    // todo:
    // get reader for fastx file
    // iterate over fastx file and store read lengths with read index
    // calculate target total length
    // initialise bitvector to length of number of reads
    // while total_length < target_total_length
    //     take rand element from lengths (using seed if given)
    //     set bit in bitvector to true
    //     add to total_length
    //     delete element from vector
    // set index of last true bit as StopIdx
    // iterate over fastx with index counter
    //     if index in bitvector is true: output read
    //     if index == StopIdx: break
    Ok(())
}
