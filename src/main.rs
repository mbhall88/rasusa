mod cli;
pub use crate::cli::Cli;
use structopt::StructOpt;

fn main() {
    let _args = Cli::from_args();
    // todo:
    // determine if file is fasta or fastq
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
}
