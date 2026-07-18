#![allow(clippy::redundant_clone)]

extern crate core;

use anyhow::Result;

pub mod alignment;
pub mod cli;
pub mod fastx;
pub mod reads;
pub mod source;
pub mod subsampler;

pub use crate::cli::Cli;
pub use crate::fastx::Fastx;
pub use crate::subsampler::SubSampler;

pub trait Runner {
    fn run(&mut self) -> Result<()>;
}
