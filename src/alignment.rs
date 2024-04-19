use std::io::stdout;
use std::path::PathBuf;
use clap::Parser;
use log::{debug, error, info, warn};
use niffler::compression;
use crate::cli::{Coverage, GenomeSize, check_path_exists, parse_fraction, parse_compression_format, parse_level, CliError};
use crate::{Cli, Fastx, Runner, SubSampler};
use anyhow::{anyhow, Context, Result};

#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct Alignment {
    /// Path to the indexed alignment file (SAM/BAM/CRAM) to subsample
    #[clap(value_parser = check_path_exists, name = "FILE")]
    pub aln: PathBuf,

    /// The desired depth of coverage to subsample the reads to
    #[clap(
    short,
    long,
    value_name = "FLOAT",
    )]
    pub coverage: Option<Coverage>,
}

impl Runner for Alignment {
    fn run(&mut self) -> Result<()> {
        info!("Subsampling alignment file: {:?}", self.aln);
        Ok(())
    }
}