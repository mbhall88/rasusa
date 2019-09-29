mod cli;

pub use crate::cli::parse_args;

fn main() {
    let args = cli::parse_args(std::env::args());
}
