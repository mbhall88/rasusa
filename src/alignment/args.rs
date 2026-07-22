use std::path::PathBuf;

use clap::{Parser, ValueEnum};
use noodles_util::alignment::io::Format;

use crate::cli::check_path_exists;
use crate::format::infer_format_from_char;

#[derive(Debug, Clone, Copy, PartialEq, ValueEnum)]
pub enum SubsamplingStrategy {
    /// A linear scan approach using sweep line algorithm with random priority. Requires sorted alignment input.
    Stream,

    /// A fetching approach to randomly subsample reads given read overlap position. Requires indexed input (.bai).
    Fetch,
}

#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct Alignment {
    /// Path to the input alignment file (SAM/BAM/CRAM) to subsample
    ///
    /// Note: An index (.bai) is required when using '--strategy fetch'.
    #[arg(value_parser = check_path_exists, name = "FILE")]
    pub aln: PathBuf,

    /// Path to the output subsampled alignment file. Defaults to stdout (same format as input)
    ///
    /// The output is not guaranteed to be sorted. We recommend piping the output to `samtools sort`
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,

    /// Output format. Rasusa will attempt to infer the format from the output file extension if not provided
    #[arg(short='O', long = "output-format", value_name = "FMT", value_parser = infer_format_from_char)]
    pub output_format: Option<Format>,

    /// The desired depth of coverage to subsample the alignment to
    #[arg(short, long, value_name = "INT", value_parser = clap::value_parser!(u32).range(1..))]
    pub coverage: u32,

    /// Random seed to use.
    #[arg(short, long, value_name = "INT")]
    pub seed: Option<u64>,

    // Strategy selection
    /// Subsampling strategy
    #[arg(long, value_enum, default_value_t = SubsamplingStrategy::Stream)]
    pub strategy: SubsamplingStrategy,

    // Algorithm specific arguments
    /// [Stream] A maximum distance (bp) allowed between start position of new read and the worst read
    /// in the heap to consider them to be 'swappable'.
    ///
    /// Larger values allow swapping reads over greater distances, but may cause local undersampling.
    /// A value of `0` means only allows swap between reads that have the same start position.
    #[arg(long, default_value_t = 5, value_name = "INT", value_parser = clap::value_parser!(i64).range(0..))]
    pub swap_distance: i64,

    /// [Fetch] When a region has less than the desired coverage, the step size to move along the chromosome
    /// to find more reads.
    ///
    /// The lowest of the step and the minimum end coordinate of the reads in the region will be used.
    /// This parameter can have a significant impact on the runtime of the subsampling process.
    #[arg(long, default_value_t = 100, value_name = "INT", value_parser = clap::value_parser!(i64).range(1..))]
    pub step_size: i64,

    /// [Fetch] The size of the genomic window (bp) to cache into memory at once.
    ///
    /// Larger values reduce disk seeking, but at the cost of high memory usage.
    /// The minimum value is 1,000 bp to avoid small region queries.
    #[arg(long, default_value_t = 10_000, value_name = "INT", value_parser = clap::value_parser!(u64).range(1000..))]
    pub batch_size: u64,
}
