use anyhow::Result;
use clap::Parser;
use env_logger::Builder;
use log::{debug, LevelFilter};

use rasusa::cli::{Cli, Commands};
use rasusa::Runner;

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
        Commands::Cite(cmd) => Box::new(cmd),
    };

    subcmd.run()?;

    Ok(())
}
