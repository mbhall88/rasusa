mod cli;

pub use crate::cli::Cli;
use structopt::StructOpt;

fn main() {
    let _args = Cli::from_args();
}
