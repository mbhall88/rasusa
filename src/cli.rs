use clap::Parser;
use regex::Regex;
use std::ffi::OsStr;
use std::fs::File;
use std::io::BufRead;
use std::ops::{Div, Mul};
use std::path::{Path, PathBuf};
use std::str::FromStr;
use thiserror::Error;

/// Randomly subsample reads to a specified coverage.
#[derive(Debug, Parser)]
#[clap(author, version, about)]
pub struct Cli {
    /// The fast{a,q} file(s) to subsample.
    ///
    /// For paired Illumina you may either pass this flag twice `-i r1.fq -i r2.fq` or give two
    /// files consecutively `-i r1.fq r2.fq`.
    #[clap(
        short = 'i',
        long = "input",
        value_parser = check_path_exists,
        num_args = 1..=2,
        required = true
    )]
    pub input: Vec<PathBuf>,

    /// Output filepath(s); stdout if not present.
    ///
    /// For paired Illumina you may either pass this flag twice `-o o1.fq -o o2.fq` or give two
    /// files consecutively `-o o1.fq o2.fq`. NOTE: The order of the pairs is assumed to be the
    /// same as that given for --input. This option is required for paired input.
    #[clap(short = 'o', long = "output", num_args = 1..=2)]
    pub output: Vec<PathBuf>,

    /// Genome size to calculate coverage with respect to. e.g., 4.3kb, 7Tb, 9000, 4.1MB
    ///
    /// Alternatively, a FASTA/Q index file can be provided and the genome size will be
    /// set to the sum of all reference sequences.
    ///
    /// If --bases is not provided, this option and --coverage are required
    #[clap(
        short,
        long,
        required_unless_present_any = &["bases", "num", "frac"],
        requires = "coverage",
        value_name = "size|faidx",
        conflicts_with_all = &["num", "frac"]
    )]
    pub genome_size: Option<GenomeSize>,

    /// The desired coverage to sub-sample the reads to
    ///
    /// If --bases is not provided, this option and --genome-size are required
    #[clap(
        short,
        long,
        value_name = "FLOAT",
        required_unless_present_any = &["bases", "num", "frac"],
        requires = "genome_size",
        conflicts_with_all = &["num", "frac"]
    )]
    pub coverage: Option<Coverage>,

    /// Explicitly set the number of bases required e.g., 4.3kb, 7Tb, 9000, 4.1MB
    ///
    /// If this option is given, --coverage and --genome-size are ignored
    #[clap(short, long, value_name = "bases", conflicts_with_all = &["num", "frac"])]
    pub bases: Option<GenomeSize>,

    /// Subsample to a specific number of reads
    ///
    /// If paired-end reads are passed, this is the number of (matched) reads from EACH file.
    /// This option accepts the same format as genome size - e.g., 1k will take 1000 reads
    #[clap(short, long, value_name = "INT", conflicts_with = "frac")]
    pub num: Option<GenomeSize>,

    /// Subsample to a fraction of the reads - e.g., 0.5 samples half the reads
    ///
    /// Values >1 and <=100 will be automatically converted - e.g., 25 => 0.25
    #[clap(short, long, value_name = "FLOAT", value_parser = parse_fraction, conflicts_with = "num")]
    pub frac: Option<f32>,

    /// Random seed to use.
    #[clap(short = 's', long = "seed", value_name = "INT")]
    pub seed: Option<u64>,

    /// Switch on verbosity.
    #[clap(short)]
    pub verbose: bool,

    /// u: uncompressed; b: Bzip2; g: Gzip; l: Lzma; x: Xz (Lzma); z: Zstd
    ///
    /// Rasusa will attempt to infer the output compression format automatically from the filename
    /// extension. This option is used to override that. If writing to stdout, the default is
    /// uncompressed
    #[clap(short = 'O', long, value_name = "u|b|g|l|x|z", value_parser = parse_compression_format)]
    pub output_type: Option<niffler::compression::Format>,

    /// Compression level to use if compressing output. Uses the default level for the format if
    /// not specified.
    #[clap(short = 'l', long, value_parser = parse_level, value_name = "1-21")]
    pub compress_level: Option<niffler::Level>,
}

impl Cli {
    /// Checks there is a valid and equal number of `--input` and `--output` arguments given.
    ///
    /// # Errors
    /// A [`CliError::BadInputOutputCombination`](#clierror) is returned for the following:
    /// - Either `--input` or `--output` are passed more than twice
    /// - An unequal number of `--input` and `--output` are passed. The only exception to
    /// this is if one `--input` and zero `--output` are passed, in which case, the output
    /// will be sent to STDOUT.
    pub fn validate_input_output_combination(&self) -> Result<(), CliError> {
        let out_len = self.output.len();
        let in_len = self.input.len();

        if in_len > 2 {
            let msg = String::from("Got more than 2 files for input.");
            return Err(CliError::BadInputOutputCombination(msg));
        }
        if out_len > 2 {
            let msg = String::from("Got more than 2 files for output.");
            return Err(CliError::BadInputOutputCombination(msg));
        }
        match in_len as isize - out_len as isize {
            diff if diff == 1 && in_len == 1 => Ok(()),
            diff if diff != 0 => Err(CliError::BadInputOutputCombination(format!(
                "Got {} --input but {} --output",
                in_len, out_len
            ))),
            _ => Ok(()),
        }
    }
}

/// A collection of custom errors relating to the command line interface for this package.
#[derive(Error, Debug, PartialEq)]
pub enum CliError {
    /// Indicates that a string cannot be parsed into a [`MetricSuffix`](#metricsuffix).
    #[error("{0} is not a valid metric suffix")]
    InvalidMetricSuffix(String),

    /// Indicates that a string cannot be parsed into a [`GenomeSize`](#genomesize).
    #[error("{0} is not a valid genome size. Valid forms include 4gb, 3000, 8.7Kb etc.")]
    InvalidGenomeSizeString(String),

    /// Indicates that a string cannot be parsed into a [`Coverage`](#coverage).
    #[error("{0} is not a valid coverage string. Coverage must be either an integer or a float and can end with an optional 'x' character")]
    InvalidCoverageValue(String),

    /// Indicates that a string cannot be parsed into a [`CompressionFormat`](#compressionformat).
    #[error("{0} is not a valid output format")]
    InvalidCompression(String),

    /// Indicates the fraction could not be parsed to the range 0-1
    #[error("{0} could not be parsed to the range 0-1")]
    FractionOutOfRange(String),

    /// Indicates a bad combination of input and output files was passed.
    #[error("Bad combination of input and output files: {0}")]
    BadInputOutputCombination(String),

    /// Faidx IO error
    #[error("Failed to open/parse faidx file {0}")]
    FaidxError(String),
}

/// A metric suffix is a unit suffix used to indicate the multiples of (in this case) base pairs.
/// For instance, the metric suffix 'Kb' refers to kilobases. Therefore, 6.9kb means 6900 base pairs.
#[derive(PartialEq, Debug)]
enum MetricSuffix {
    Base,
    Kilo,
    Mega,
    Giga,
    Tera,
}

impl FromStr for MetricSuffix {
    type Err = CliError;

    /// Parses a string into a `MetricSuffix`.
    ///
    /// # Example
    /// ```rust
    /// let s = "5.5mb";
    /// let metric_suffix = MetricSuffix::from_str(s);
    ///
    /// assert_eq!(metric_suffix, MetricSuffix::Mega)
    /// ```
    fn from_str(suffix: &str) -> Result<Self, Self::Err> {
        let suffix_lwr = suffix.to_lowercase();
        let metric_suffix = match suffix_lwr.as_str() {
            s if "b".contains(s) => MetricSuffix::Base,
            s if "kb".contains(s) => MetricSuffix::Kilo,
            s if "mb".contains(s) => MetricSuffix::Mega,
            s if "gb".contains(s) => MetricSuffix::Giga,
            s if "tb".contains(s) => MetricSuffix::Tera,
            _ => {
                return Err(CliError::InvalidMetricSuffix(suffix.to_string()));
            }
        };
        Ok(metric_suffix)
    }
}

/// Allow for multiplying a `f64` by a `MetricSuffix`.
///
/// # Example
///
/// ```rust
/// let metric_suffix = MetricSuffix::Mega;
/// let x: f64 = 5.5;
///
/// assert_eq!(x * metric_suffix, 5_500_000)
/// ```
impl Mul<MetricSuffix> for f64 {
    type Output = Self;

    fn mul(self, rhs: MetricSuffix) -> Self::Output {
        match rhs {
            MetricSuffix::Base => self,
            MetricSuffix::Kilo => self * 1_000.0,
            MetricSuffix::Mega => self * 1_000_000.0,
            MetricSuffix::Giga => self * 1_000_000_000.0,
            MetricSuffix::Tera => self * 1_000_000_000_000.0,
        }
    }
}

/// An object for collecting together methods for working with the genome size parameter for this
/// package.
#[derive(Debug, PartialOrd, PartialEq, Copy, Clone)]
pub struct GenomeSize(u64);

/// Allow for comparison of a `u64` and a `GenomeSize`.
///
/// # Example
///
/// ```rust
/// assert!(GenomeSize(10) == 10)
/// ```
impl PartialEq<u64> for GenomeSize {
    fn eq(&self, other: &u64) -> bool {
        self.0 == *other
    }
}

impl From<GenomeSize> for u64 {
    fn from(g: GenomeSize) -> Self {
        g.0
    }
}

// add method to convert to f64
impl From<GenomeSize> for f64 {
    fn from(g: GenomeSize) -> Self {
        g.0 as f64
    }
}

impl FromStr for GenomeSize {
    type Err = CliError;

    /// Parses a string into a `GenomeSize`.
    ///
    /// # Example
    /// ```rust
    /// let s = "5.5mb";
    /// let genome_size = GenomeSize::from_str(s);
    ///
    /// assert_eq!(genome_size, GenomeSize(5_500_000))
    /// ```
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let text = s.to_lowercase();
        let re = Regex::new(r"(?P<size>[0-9]*\.?[0-9]+)(?P<sfx>\w*)$").unwrap();
        let captures = match re.captures(text.as_str()) {
            Some(cap) => cap,
            None => match Self::from_faidx(Path::new(&s)) {
                Ok(g) => return Ok(g),
                _ => return Err(CliError::InvalidGenomeSizeString(s.to_string())),
            },
        };
        let size = captures
            .name("size")
            .unwrap()
            .as_str()
            .parse::<f64>()
            .unwrap();
        let metric_suffix = MetricSuffix::from_str(captures.name("sfx").unwrap().as_str())?;

        Ok(GenomeSize((size * metric_suffix) as u64))
    }
}

impl GenomeSize {
    fn from_faidx(path: &Path) -> Result<Self, CliError> {
        let mut size = 0;
        let file = File::open(path)
            .map_err(|_| CliError::FaidxError(path.to_string_lossy().to_string()))?;
        let rdr = std::io::BufReader::new(file);
        for result in rdr.lines() {
            let line =
                result.map_err(|_| CliError::FaidxError(path.to_string_lossy().to_string()))?;
            let fields: Vec<&str> = line.split('\t').collect();
            size += u64::from_str(fields[1])
                .map_err(|_| CliError::FaidxError(path.to_string_lossy().to_string()))?;
        }
        Ok(GenomeSize(size))
    }
}

/// Allow for multiplying a `GenomeSize` by a [`Coverage`](#coverage).
///
/// # Example
///
/// ```rust
/// let genome_size = GenomeSize(100);
/// let covg = Coverage(5);
///
/// assert_eq!(genome_size * covg, 500)
/// ```
impl Mul<Coverage> for GenomeSize {
    type Output = u64;

    fn mul(self, rhs: Coverage) -> Self::Output {
        (self.0 as f32 * rhs.0) as u64
    }
}

/// Allow for dividing a `u64` by a `GenomeSize`.
///
/// # Example
///
/// ```rust
/// let x: u64 = 210;
/// let size = GenomeSize(200);
///
/// let actual = x / size;
/// let expected = 1.05;
///
/// assert_eq!(actual, expected)
/// ```
impl Div<GenomeSize> for u64 {
    type Output = f64;

    fn div(self, rhs: GenomeSize) -> Self::Output {
        (self as f64) / (rhs.0 as f64)
    }
}

/// An object for collecting together methods for working with the coverage parameter for this
/// package.
#[derive(Debug, PartialOrd, PartialEq, Copy, Clone)]
pub struct Coverage(pub f32);

/// Allow for comparison of a `f32` and a `Coverage`.
///
/// # Example
///
/// ```rust
/// assert!(Coverage(10) == 10.0)
/// ```
impl PartialEq<f32> for Coverage {
    fn eq(&self, other: &f32) -> bool {
        self.0 == *other
    }
}

impl FromStr for Coverage {
    type Err = CliError;

    /// Parses a string into a `Coverage`.
    ///
    /// # Example
    /// ```rust
    /// let s = "100x";
    /// let covg = Coverage::from_str(s);
    ///
    /// assert_eq!(covg, Coverage(100))
    /// ```
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let re = Regex::new(r"^(?P<covg>[0-9]*\.?[0-9]+)(?i)x?$").unwrap();
        let captures = match re.captures(s) {
            Some(cap) => cap,
            None => {
                return Err(CliError::InvalidCoverageValue(s.to_string()));
            }
        };
        Ok(Coverage(
            captures
                .name("covg")
                .unwrap()
                .as_str()
                .parse::<f32>()
                .unwrap(),
        ))
    }
}

/// Allow for multiplying a `Coverage` by a [`GenomeSize`](#genomesize).
///
/// # Example
///
/// ```rust
/// let covg = Coverage(5);
/// let genome_size = GenomeSize(100);
///
/// assert_eq!(covg * genome_size, 500)
/// ```
impl Mul<GenomeSize> for Coverage {
    type Output = u64;

    fn mul(self, rhs: GenomeSize) -> Self::Output {
        (self.0 * (rhs.0 as f32)) as u64
    }
}

pub trait CompressionExt {
    fn from_path<S: AsRef<OsStr> + ?Sized>(p: &S) -> Self;
}

impl CompressionExt for niffler::compression::Format {
    /// Attempts to infer the compression type from the file extension. If the extension is not
    /// known, then Uncompressed is returned.
    fn from_path<S: AsRef<OsStr> + ?Sized>(p: &S) -> Self {
        let path = Path::new(p);
        match path.extension().map(|s| s.to_str()) {
            Some(Some("gz")) => Self::Gzip,
            Some(Some("bz") | Some("bz2")) => Self::Bzip,
            Some(Some("lzma")) => Self::Lzma,
            _ => Self::No,
        }
    }
}

fn parse_compression_format(s: &str) -> Result<niffler::compression::Format, CliError> {
    match s {
        "b" | "B" => Ok(niffler::Format::Bzip),
        "g" | "G" => Ok(niffler::Format::Gzip),
        "l" | "L" => Ok(niffler::Format::Lzma),
        "x" | "X" => Ok(niffler::Format::Xz),
        "z" | "Z" => Ok(niffler::Format::Zstd),
        "u" | "U" => Ok(niffler::Format::No),
        _ => Err(CliError::InvalidCompression(s.to_string())),
    }
}
/// A utility function that allows the CLI to error if a path doesn't exist
fn check_path_exists<S: AsRef<OsStr> + ?Sized>(s: &S) -> Result<PathBuf, String> {
    let path = PathBuf::from(s);
    if path.exists() {
        Ok(path)
    } else {
        Err(format!("{:?} does not exist", path))
    }
}

/// A utility function to validate compression level is in allowed range
#[allow(clippy::redundant_clone)]
fn parse_level(s: &str) -> Result<niffler::Level, String> {
    let lvl = match s.parse::<u8>() {
        Ok(1) => niffler::Level::One,
        Ok(2) => niffler::Level::Two,
        Ok(3) => niffler::Level::Three,
        Ok(4) => niffler::Level::Four,
        Ok(5) => niffler::Level::Five,
        Ok(6) => niffler::Level::Six,
        Ok(7) => niffler::Level::Seven,
        Ok(8) => niffler::Level::Eight,
        Ok(9) => niffler::Level::Nine,
        Ok(10) => niffler::Level::Ten,
        Ok(11) => niffler::Level::Eleven,
        Ok(12) => niffler::Level::Twelve,
        Ok(13) => niffler::Level::Thirteen,
        Ok(14) => niffler::Level::Fourteen,
        Ok(15) => niffler::Level::Fifteen,
        Ok(16) => niffler::Level::Sixteen,
        Ok(17) => niffler::Level::Seventeen,
        Ok(18) => niffler::Level::Eighteen,
        Ok(19) => niffler::Level::Nineteen,
        Ok(20) => niffler::Level::Twenty,
        Ok(21) => niffler::Level::TwentyOne,
        _ => return Err(format!("Compression level {} not in the range 1-21", s)),
    };
    Ok(lvl)
}

/// A utility function to parse the fraction CLI option to the range 0-1
fn parse_fraction(s: &str) -> Result<f32, CliError> {
    let f = f32::from_str(s).map_err(|_| CliError::FractionOutOfRange(s.to_string()))?;
    if !(0.0..=100.0).contains(&f) {
        return Err(CliError::FractionOutOfRange(s.to_string()));
    }

    let result = if f.ceil() > 1.0_f32 { f / 100.0_f32 } else { f };

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    const ERROR: f32 = f32::EPSILON;

    #[test]
    fn parse_fraction_negative_fraction_returns_error() {
        let s = "-0.1";

        let actual = parse_fraction(s).unwrap_err();
        let expected = CliError::FractionOutOfRange(s.to_string());

        assert_eq!(actual, expected)
    }

    #[test]
    fn parse_fraction_zero_is_zero() {
        let s = "0";

        let actual = parse_fraction(s).unwrap();
        let expected = 0.0_f32;

        assert!((expected - actual).abs() < ERROR)
    }

    #[test]
    fn parse_fraction_one_is_one() {
        let s = "1";

        let actual = parse_fraction(s).unwrap();
        let expected = 1.0_f32;

        assert!((expected - actual).abs() < ERROR)
    }

    #[test]
    fn parse_fraction_one_asf32_is_one() {
        let s = "1.0";

        let actual = parse_fraction(s).unwrap();
        let expected = 1.0_f32;

        assert!((expected - actual).abs() < ERROR)
    }

    #[test]
    fn parse_fraction_less_one_unchanged() {
        let s = "0.25";

        let actual = parse_fraction(s).unwrap();
        let expected = 0.25_f32;

        assert!((expected - actual).abs() < ERROR)
    }

    #[test]
    fn parse_fraction_over_100_fails() {
        let s = "101.5";

        let actual = parse_fraction(s).unwrap_err();
        let expected = CliError::FractionOutOfRange(s.to_string());

        assert_eq!(actual, expected)
    }

    #[test]
    fn parse_fraction_check_upper_range() {
        let s = "100";

        let actual = parse_fraction(s).unwrap();
        let expected = 1.0_f32;

        assert!((expected - actual).abs() < ERROR)
    }

    #[test]
    fn parse_fraction_check_lower_range() {
        let s = "1";

        let actual = parse_fraction(s).unwrap();
        let expected = 1_f32;

        assert!((expected - actual).abs() < ERROR)
    }

    #[test]
    fn parse_fraction_non_fraction_fails() {
        let s = "foo";

        let actual = parse_fraction(s).unwrap_err();
        let expected = CliError::FractionOutOfRange(s.to_string());

        assert_eq!(actual, expected)
    }

    #[test]
    fn parse_fraction_over_one_converts() {
        let s = "25";

        let actual = parse_fraction(s).unwrap();
        let expected = 0.25_f32;

        assert!((expected - actual).abs() < ERROR)
    }

    #[test]
    fn parse_fraction_over_one_as_fration_converts() {
        let s = "25.6";

        let actual = parse_fraction(s).unwrap();
        let expected = 0.256_f32;

        assert!((expected - actual).abs() < ERROR)
    }

    #[test]
    fn no_args_given_raises_error() {
        let passed_args = vec!["rasusa"];
        let args: Result<Cli, clap::Error> = Cli::try_parse_from(passed_args);

        let actual = args.unwrap_err().kind();
        let expected = clap::error::ErrorKind::MissingRequiredArgument;

        assert_eq!(actual, expected)
    }

    #[test]
    fn no_input_file_given_raises_error() {
        let passed_args = vec!["rasusa", "-c", "30", "-g", "3mb"];
        let args: Result<Cli, clap::Error> = Cli::try_parse_from(passed_args);

        let actual = args.unwrap_err().kind();
        let expected = clap::error::ErrorKind::MissingRequiredArgument;

        assert_eq!(actual, expected)
    }

    #[test]
    fn no_coverage_given_raises_error() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec!["rasusa", "-i", infile, "-g", "3mb"];
        let args: Result<Cli, clap::Error> = Cli::try_parse_from(passed_args);

        let actual = args.unwrap_err().kind();
        let expected = clap::error::ErrorKind::MissingRequiredArgument;

        assert_eq!(actual, expected)
    }

    #[test]
    fn invalid_coverage_given_raises_error() {
        let passed_args = vec!["rasusa", "-i", "in.fq", "-g", "3mb", "-c", "foo"];
        let args: Result<Cli, clap::Error> = Cli::try_parse_from(passed_args);

        let actual = args.unwrap_err().kind();
        let expected = clap::error::ErrorKind::ValueValidation;

        assert_eq!(actual, expected)
    }

    #[test]
    fn no_genome_size_given_raises_error() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec!["rasusa", "-i", infile, "-c", "5"];
        let args: Result<Cli, clap::Error> = Cli::try_parse_from(passed_args);

        let actual = args.unwrap_err().kind();
        let expected = clap::error::ErrorKind::MissingRequiredArgument;

        assert_eq!(actual, expected)
    }

    #[test]
    fn no_genome_size_but_bases_and_coverage() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec!["rasusa", "-i", infile, "-b", "5", "-c", "7"];

        let args: Result<Cli, clap::Error> = Cli::try_parse_from(passed_args);

        let actual = args.unwrap_err().kind();
        let expected = clap::error::ErrorKind::MissingRequiredArgument;

        assert_eq!(actual, expected)
    }

    #[test]
    fn bases_and_coverage_and_genome_size_all_allowed() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec!["rasusa", "-i", infile, "-b", "5", "-c", "7", "-g", "5m"];

        let args = Cli::try_parse_from(passed_args).unwrap();

        assert_eq!(args.bases.unwrap(), GenomeSize(5));
    }

    #[test]
    fn no_genome_size_or_coverage_given_but_bases_given() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec!["rasusa", "-i", infile, "-b", "5"];

        let args = Cli::try_parse_from(passed_args).unwrap();

        assert_eq!(args.bases.unwrap(), GenomeSize(5));
    }

    #[test]
    fn no_genome_size_or_coverage_given_but_num_given() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec!["rasusa", "-i", infile, "-n", "5"];

        let args = Cli::try_parse_from(passed_args).unwrap();

        assert_eq!(args.num.unwrap(), GenomeSize(5));
    }

    #[test]
    fn no_genome_size_or_coverage_given_but_frac_given() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec!["rasusa", "-i", infile, "-f", "5"];

        let args = Cli::try_parse_from(passed_args).unwrap();

        assert!((args.frac.unwrap() - 0.05_f32).abs() < ERROR);
    }

    #[test]
    fn num_and_coverage_and_genome_size_not_allowed() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec!["rasusa", "-i", infile, "-n", "5", "-c", "7", "-g", "5m"];

        let args: Result<Cli, clap::Error> = Cli::try_parse_from(passed_args);

        let actual = args.unwrap_err().kind();
        let expected = clap::error::ErrorKind::ArgumentConflict;

        assert_eq!(actual, expected)
    }

    #[test]
    fn frac_and_coverage_and_genome_size_not_allowed() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec!["rasusa", "-i", infile, "-f", "5", "-c", "7", "-g", "5m"];

        let args: Result<Cli, clap::Error> = Cli::try_parse_from(passed_args);

        let actual = args.unwrap_err().kind();
        let expected = clap::error::ErrorKind::ArgumentConflict;

        assert_eq!(actual, expected)
    }

    #[test]
    fn frac_and_bases_not_allowed() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec!["rasusa", "-i", infile, "-f", "5", "-b", "7"];

        let args: Result<Cli, clap::Error> = Cli::try_parse_from(passed_args);

        let actual = args.unwrap_err().kind();
        let expected = clap::error::ErrorKind::ArgumentConflict;

        assert_eq!(actual, expected)
    }

    #[test]
    fn frac_and_num_not_allowed() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec!["rasusa", "-i", infile, "-f", "5", "-n", "7"];

        let args: Result<Cli, clap::Error> = Cli::try_parse_from(passed_args);

        let actual = args.unwrap_err().kind();
        let expected = clap::error::ErrorKind::ArgumentConflict;

        assert_eq!(actual, expected)
    }

    #[test]
    fn bases_and_num_not_allowed() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec!["rasusa", "-i", infile, "-b", "5", "-n", "7"];

        let args: Result<Cli, clap::Error> = Cli::try_parse_from(passed_args);

        let actual = args.unwrap_err().kind();
        let expected = clap::error::ErrorKind::ArgumentConflict;

        assert_eq!(actual, expected)
    }

    #[test]
    fn invalid_genome_size_given_raises_error() {
        let passed_args = vec!["rasusa", "-i", "in.fq", "-c", "5", "-g", "8jb"];
        let args: Result<Cli, clap::Error> = Cli::try_parse_from(passed_args);

        let actual = args.unwrap_err().kind();
        let expected = clap::error::ErrorKind::ValueValidation;

        assert_eq!(actual, expected)
    }

    #[test]
    fn faidx_given_as_genome_size() {
        let infile = "tests/cases/r1.fq.gz";
        let faidx = "tests/cases/h37rv.fa.fai";
        let passed_args = vec!["rasusa", "-i", infile, "-c", "5", "-g", faidx];
        let args = Cli::try_parse_from(passed_args).unwrap();

        assert_eq!(args.genome_size.unwrap(), GenomeSize(4411532))
    }

    #[test]
    fn invalid_seed_given_raises_error() {
        let passed_args = vec!["rasusa", "-i", "in.fq", "-c", "5", "-g", "8mb", "-s", "foo"];
        let args: Result<Cli, clap::Error> = Cli::try_parse_from(passed_args);

        let actual = args.unwrap_err().kind();
        let expected = clap::error::ErrorKind::ValueValidation;

        assert_eq!(actual, expected)
    }

    #[test]
    fn all_valid_args_parsed_as_expected() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![
            "rasusa",
            "-i",
            infile,
            "-c",
            "5",
            "-g",
            "8mb",
            "-s",
            "88",
            "-o",
            "my/output/file.fq",
        ];
        let args = Cli::try_parse_from(passed_args).unwrap();

        assert_eq!(args.input[0], PathBuf::from_str(infile).unwrap());
        assert_eq!(args.coverage.unwrap(), Coverage(5.0));
        assert_eq!(args.genome_size.unwrap(), GenomeSize(8_000_000));
        assert_eq!(args.seed, Some(88));
        assert_eq!(
            args.output[0],
            PathBuf::from_str("my/output/file.fq").unwrap()
        )
    }

    #[test]
    fn all_valid_args_with_two_inputs_parsed_as_expected() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![
            "rasusa",
            "-i",
            infile,
            infile,
            "-c",
            "5",
            "-g",
            "8mb",
            "-s",
            "88",
            "-o",
            "my/output/file.fq",
        ];
        let args = Cli::try_parse_from(passed_args).unwrap();

        let expected_input = vec![
            PathBuf::from_str(infile).unwrap(),
            PathBuf::from_str(infile).unwrap(),
        ];
        assert_eq!(args.input, expected_input);
        assert_eq!(args.coverage.unwrap(), Coverage(5.0));
        assert_eq!(args.genome_size.unwrap(), GenomeSize(8_000_000));
        assert_eq!(args.seed, Some(88));
        assert_eq!(
            args.output[0],
            PathBuf::from_str("my/output/file.fq").unwrap()
        )
    }

    #[test]
    fn all_valid_args_with_two_inputs_using_flag_twice_parsed_as_expected() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![
            "rasusa",
            "-i",
            infile,
            "-i",
            infile,
            "-c",
            "5",
            "-g",
            "8mb",
            "-s",
            "88",
            "-o",
            "my/output/file.fq",
        ];
        let args = Cli::try_parse_from(passed_args).unwrap();

        let expected_input = vec![
            PathBuf::from_str(infile).unwrap(),
            PathBuf::from_str(infile).unwrap(),
        ];
        assert_eq!(args.input, expected_input);
        assert_eq!(args.coverage.unwrap(), Coverage(5.0));
        assert_eq!(args.genome_size.unwrap(), GenomeSize(8_000_000));
        assert_eq!(args.seed, Some(88));
        assert_eq!(
            args.output[0],
            PathBuf::from_str("my/output/file.fq").unwrap()
        )
    }

    #[test]
    fn three_inputs_raises_error() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![
            "rasusa",
            "-i",
            infile,
            "-i",
            infile,
            "-i",
            infile,
            "-c",
            "5",
            "-g",
            "8mb",
            "-s",
            "88",
            "-o",
            "my/output/file.fq",
        ];
        let args = Cli::try_parse_from(passed_args).unwrap();

        let actual: CliError = args.validate_input_output_combination().unwrap_err();
        let expected =
            CliError::BadInputOutputCombination(String::from("Got more than 2 files for input."));

        assert_eq!(actual, expected)
    }

    #[test]
    fn three_outputs_raises_error() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![
            "rasusa",
            "-i",
            infile,
            "-i",
            infile,
            "-c",
            "5",
            "-g",
            "8mb",
            "-s",
            "88",
            "-o",
            "my/output/file.fq",
            "-o",
            "out.fq",
            "-o",
            "out.fq",
        ];
        let args = Cli::try_parse_from(passed_args).unwrap();

        let actual: CliError = args.validate_input_output_combination().unwrap_err();
        let expected =
            CliError::BadInputOutputCombination(String::from("Got more than 2 files for output."));

        assert_eq!(actual, expected)
    }

    #[test]
    fn one_input_no_output_is_ok() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec!["rasusa", "-i", infile, "-c", "5", "-g", "8mb", "-s", "88"];
        let args = Cli::try_parse_from(passed_args).unwrap();

        let actual = args.validate_input_output_combination();

        assert!(actual.is_ok())
    }

    #[test]
    fn two_inputs_one_output_raises_error() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![
            "rasusa", "-i", infile, "-i", infile, "-c", "5", "-g", "8mb", "-s", "88", "-o",
            "out.fq",
        ];
        let args = Cli::try_parse_from(passed_args).unwrap();

        let actual: CliError = args.validate_input_output_combination().unwrap_err();
        let expected =
            CliError::BadInputOutputCombination(String::from("Got 2 --input but 1 --output"));

        assert_eq!(actual, expected)
    }

    #[test]
    fn one_input_two_outputs_raises_error() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![
            "rasusa", "-i", infile, "-c", "5", "-g", "8mb", "-s", "88", "-o", "out.fq", "-o",
            "out.fq",
        ];
        let args = Cli::try_parse_from(passed_args).unwrap();

        let actual: CliError = args.validate_input_output_combination().unwrap_err();
        let expected =
            CliError::BadInputOutputCombination(String::from("Got 1 --input but 2 --output"));

        assert_eq!(actual, expected)
    }

    #[test]
    fn two_input_two_outputs_is_ok() {
        let infile = "tests/cases/r1.fq.gz";
        let passed_args = vec![
            "rasusa", "-i", infile, "-i", infile, "-c", "5", "-g", "8mb", "-s", "88", "-o",
            "out.fq", "-o", "out.fq",
        ];
        let args = Cli::try_parse_from(passed_args).unwrap();

        let actual = args.validate_input_output_combination();

        assert!(actual.is_ok())
    }

    #[test]
    fn float_multiply_with_base_unit() {
        let actual = 4.5 * MetricSuffix::Base;
        let expected = 4.5;

        let diff = (actual.abs() - expected).abs();
        assert!(diff < f64::EPSILON)
    }

    #[test]
    fn integer_only_returns_integer() {
        let actual = GenomeSize::from_str("6").unwrap();
        let expected = 6;

        assert_eq!(actual, expected);
    }

    #[test]
    fn float_only_returns_integer() {
        let actual = GenomeSize::from_str("6.5").unwrap();
        let expected = 6;

        assert_eq!(actual, expected);
    }

    #[test]
    fn int_and_suffix_returns_multiplied_int() {
        let actual = GenomeSize::from_str("5mb").unwrap();
        let expected = 5_000_000;

        assert_eq!(actual, expected);
    }

    #[test]
    fn float_and_suffix_returns_multiplied_float_as_int() {
        let actual = GenomeSize::from_str("5.4kB").unwrap();
        let expected = 5_400;

        assert_eq!(actual, expected);
    }

    #[test]
    fn float_without_leading_int_and_suffix_returns_multiplied_float_as_int() {
        let actual = GenomeSize::from_str(".77G").unwrap();
        let expected = 770_000_000;

        assert_eq!(actual, expected);
    }

    #[test]
    fn int_and_tera_suffix_returns_multiplied_int() {
        let actual = GenomeSize::from_str("7TB").unwrap();
        let expected = 7_000_000_000_000;

        assert_eq!(actual, expected);
    }

    #[test]
    fn int_and_base_suffix_returns_int_without_scaling() {
        let actual = GenomeSize::from_str("7B").unwrap();
        let expected = 7;

        assert_eq!(actual, expected);
    }

    #[test]
    fn invalid_suffix_returns_err() {
        let genome_size = String::from(".77uB");
        let actual = GenomeSize::from_str(genome_size.as_str()).unwrap_err();
        let expected = CliError::InvalidMetricSuffix(String::from("ub"));

        assert_eq!(actual, expected);
    }

    #[test]
    fn empty_string_returns_error() {
        let actual = GenomeSize::from_str("").unwrap_err();
        let expected = CliError::InvalidGenomeSizeString(String::from(""));

        assert_eq!(actual, expected);
    }

    #[test]
    fn suffix_with_no_size_returns_error() {
        let actual = GenomeSize::from_str("gb");

        assert!(actual.is_err());
    }

    #[test]
    fn genome_size_from_faidx() {
        let p = Path::new("tests/cases/h37rv.fa.fai");
        let actual = GenomeSize::from_faidx(p).unwrap();
        let expected = GenomeSize(4411532);

        assert_eq!(actual, expected)
    }

    #[test]
    fn genome_size_from_faidx_fastq_index() {
        let p = Path::new("tests/cases/file1.fq.fai");
        let actual = GenomeSize::from_faidx(p).unwrap();
        let expected = GenomeSize(10050);

        assert_eq!(actual, expected)
    }

    #[test]
    fn genome_size_from_str_fastq_index() {
        let actual = GenomeSize::from_str("tests/cases/file1.fq.fai").unwrap();
        let expected = GenomeSize(10050);

        assert_eq!(actual, expected)
    }

    #[test]
    fn int_coverage_returns_float() {
        let actual = Coverage::from_str("56").unwrap();
        let expected = 56.0;

        assert!((expected - actual.0).abs() < ERROR)
    }

    #[test]
    fn float_coverage_returns_float() {
        let actual = Coverage::from_str("56.6").unwrap();
        let expected = 56.6;

        assert!((expected - actual.0).abs() < ERROR)
    }

    #[test]
    fn empty_coverage_returns_err() {
        let coverage = String::from("");

        let actual = Coverage::from_str(coverage.as_str()).unwrap_err();
        let expected = CliError::InvalidCoverageValue(coverage);

        assert_eq!(actual, expected)
    }

    #[test]
    fn non_number_coverage_returns_err() {
        let coverage = String::from("foo");

        let actual = Coverage::from_str(coverage.as_str()).unwrap_err();
        let expected = CliError::InvalidCoverageValue(coverage);

        assert_eq!(actual, expected)
    }

    #[test]
    fn zero_coverage_returns_zero() {
        let actual = Coverage::from_str("0").unwrap();
        let expected = 0.0;

        assert!((expected - actual.0).abs() < ERROR)
    }

    #[test]
    fn int_ending_in_x_coverage_returns_float() {
        let actual = Coverage::from_str("1X").unwrap();
        let expected = 1.0;

        assert!((expected - actual.0).abs() < ERROR)
    }

    #[test]
    fn float_ending_in_x_coverage_returns_float() {
        let actual = Coverage::from_str("1.9X").unwrap();
        let expected = 1.9;

        assert!((expected - actual.0).abs() < ERROR)
    }

    #[test]
    fn mega_suffix_from_string() {
        let actual = MetricSuffix::from_str("MB").unwrap();
        let expected = MetricSuffix::Mega;

        assert_eq!(actual, expected)
    }

    #[test]
    fn kilo_suffix_from_string() {
        let actual = MetricSuffix::from_str("kB").unwrap();
        let expected = MetricSuffix::Kilo;

        assert_eq!(actual, expected)
    }

    #[test]
    fn giga_suffix_from_string() {
        let actual = MetricSuffix::from_str("Gb").unwrap();
        let expected = MetricSuffix::Giga;

        assert_eq!(actual, expected)
    }

    #[test]
    fn tera_suffix_from_string() {
        let actual = MetricSuffix::from_str("tb").unwrap();
        let expected = MetricSuffix::Tera;

        assert_eq!(actual, expected)
    }

    #[test]
    fn base_suffix_from_string() {
        let actual = MetricSuffix::from_str("B").unwrap();
        let expected = MetricSuffix::Base;

        assert_eq!(actual, expected)
    }

    #[test]
    fn empty_string_is_base_metric_suffix() {
        let suffix = String::from("");
        let actual = MetricSuffix::from_str(suffix.as_str()).unwrap();
        let expected = MetricSuffix::Base;

        assert_eq!(actual, expected)
    }

    #[test]
    fn invalid_suffix_raises_error() {
        let suffix = String::from("ub");
        let actual = MetricSuffix::from_str(suffix.as_str()).unwrap_err();
        let expected = CliError::InvalidMetricSuffix(suffix);

        assert_eq!(actual, expected)
    }

    #[test]
    fn multiply_genome_size_by_coverage() {
        let genome_size = GenomeSize::from_str("4.2kb").unwrap();
        let covg = Coverage::from_str("11.7866").unwrap();

        let actual = genome_size * covg;
        let expected: u64 = 49_503;

        assert_eq!(actual, expected)
    }

    #[test]
    fn multiply_coverage_by_genome_size() {
        let genome_size = GenomeSize::from_str("4.2kb").unwrap();
        let covg = Coverage::from_str("11.7866").unwrap();

        let actual = covg * genome_size;
        let expected: u64 = 49_503;

        assert_eq!(actual, expected)
    }

    #[test]
    fn divide_u64_by_genome_size() {
        let x: u64 = 210;
        let size = GenomeSize(200);

        let actual = x / size;
        let expected = 1.05;

        let diff = (actual.abs() - expected).abs();
        assert!(diff < f64::EPSILON)
    }

    #[test]
    fn compression_format_from_str() {
        let mut s = "B";
        assert_eq!(parse_compression_format(s).unwrap(), niffler::Format::Bzip);

        s = "g";
        assert_eq!(parse_compression_format(s).unwrap(), niffler::Format::Gzip);

        s = "l";
        assert_eq!(parse_compression_format(s).unwrap(), niffler::Format::Lzma);

        s = "U";
        assert_eq!(parse_compression_format(s).unwrap(), niffler::Format::No);

        s = "a";
        assert_eq!(
            parse_compression_format(s).unwrap_err(),
            CliError::InvalidCompression(s.to_string())
        );
    }

    #[test]
    fn test_in_compress_range() {
        assert!(parse_level("1").is_ok());
        assert!(parse_level("9").is_ok());
        assert!(parse_level("0").is_err());
        assert!(parse_level("22").is_err());
        assert!(parse_level("f").is_err());
        assert!(parse_level("5.5").is_err());
        assert!(parse_level("-3").is_err());
    }

    #[test]
    fn compression_format_from_path() {
        assert_eq!(niffler::Format::from_path("foo.gz"), niffler::Format::Gzip);
        assert_eq!(
            niffler::Format::from_path(Path::new("foo.gz")),
            niffler::Format::Gzip
        );
        assert_eq!(niffler::Format::from_path("baz"), niffler::Format::No);
        assert_eq!(niffler::Format::from_path("baz.fq"), niffler::Format::No);
        assert_eq!(
            niffler::Format::from_path("baz.fq.bz2"),
            niffler::Format::Bzip
        );
        assert_eq!(
            niffler::Format::from_path("baz.fq.bz"),
            niffler::Format::Bzip
        );
        assert_eq!(
            niffler::Format::from_path("baz.fq.lzma"),
            niffler::Format::Lzma
        );
    }
}
