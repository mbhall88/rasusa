use crate::cli::{check_path_exists, Coverage};
use crate::Runner;
use anyhow::Result;
use clap::Parser;
use log::info;
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct Alignment {
    /// Path to the indexed alignment file (SAM/BAM/CRAM) to subsample
    #[clap(value_parser = check_path_exists, name = "FILE")]
    pub aln: PathBuf,

    /// The desired depth of coverage to subsample the reads to
    #[clap(short, long, value_name = "FLOAT")]
    pub coverage: Option<Coverage>,
}

impl Runner for Alignment {
    fn run(&mut self) -> Result<()> {
        info!("Subsampling alignment file: {:?}", self.aln);
        Ok(())
    }
}
