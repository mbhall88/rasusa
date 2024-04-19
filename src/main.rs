#![allow(clippy::redundant_clone)]

extern crate core;

// due to a structopt problem
use std::io::stdout;

use anyhow::{Context, Result};
use clap::Parser;
use env_logger::Builder;
use log::{debug, error, info, warn, LevelFilter};
use niffler::compression;

pub use crate::cli::Cli;
use crate::cli::{Commands, Coverage};
pub use crate::fastx::Fastx;
pub use crate::subsampler::SubSampler;

mod cli;
mod fastx;
mod subsampler;
mod reads;
mod alignment;

pub trait Runner {
    fn run(&mut self) -> Result<()>;
}

fn main() -> Result<()> {
    let args: Cli = Cli::parse();
    // Initialize logger
    let log_lvl = if args.verbose {
        LevelFilter::Debug
    } else {
        LevelFilter::Info
    };
    let mut log_builder = Builder::new();
    log_builder
        .filter(None, log_lvl)
        .format_module_path(false)
        .format_target(false)
        .init();

    debug!("{:?}", args);

    let mut subcmd: Box<dyn Runner> = match args.command {
        Commands::Reads(cmd) => Box::new(cmd),
        Commands::Alignment(cmd) => Box::new(cmd),
    };

    subcmd.run()?;

    Ok(())
}
