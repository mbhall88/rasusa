use crate::alignment::Alignment;
use crate::reads::Reads;
use crate::Runner;
use clap::{Parser, Subcommand};
use regex::Regex;
use std::ffi::OsStr;
use std::fs::File;
use std::io::BufRead;
use std::ops::{Div, Mul};
use std::path::{Path, PathBuf};
use std::str::FromStr;
use thiserror::Error;

const CITATION: &str = r#"@article{Hall2022,
  doi = {10.21105/joss.03941},
  url = {https://doi.org/10.21105/joss.03941},
  year = {2022},
  publisher = {The Open Journal},
  volume = {7},
  number = {69},
  pages = {3941},
  author = {Michael B. Hall},
  title = {Rasusa: Randomly subsample sequencing reads to a specified coverage},
  journal = {Journal of Open Source Software}
}"#;

/// Randomly subsample reads or alignments
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(propagate_version = true)]
pub struct Cli {
    /// Switch on verbosity.
    #[arg(short)]
    pub verbose: bool,

    #[command(subcommand)]
    pub(crate) command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Randomly subsample reads
    Reads(Reads),
    /// Randomly subsample alignments to a specified depth of coverage
    #[command(name = "aln")]
    Alignment(Alignment),
    /// Get a bibtex formatted citation for this package.
    Cite(Cite),
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

#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct Cite {}

impl Runner for Cite {
    fn run(&mut self) -> anyhow::Result<()> {
        println!("{}", CITATION);
        Ok(())
    }
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

pub(crate) fn parse_compression_format(s: &str) -> Result<niffler::compression::Format, CliError> {
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
pub(crate) fn check_path_exists<S: AsRef<OsStr> + ?Sized>(s: &S) -> Result<PathBuf, String> {
    let path = PathBuf::from(s);
    if path.exists() {
        Ok(path)
    } else {
        Err(format!("{:?} does not exist", path))
    }
}

/// A utility function to validate compression level is in allowed range
#[allow(clippy::redundant_clone)]
pub(crate) fn parse_level(s: &str) -> Result<niffler::Level, String> {
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
pub(crate) fn parse_fraction(s: &str) -> Result<f32, CliError> {
    let f = f32::from_str(s).map_err(|_| CliError::FractionOutOfRange(s.to_string()))?;
    if !(0.0..=100.0).contains(&f) {
        return Err(CliError::FractionOutOfRange(s.to_string()));
    }

    let result = if f.ceil() > 1.0_f32 { f / 100.0_f32 } else { f };

    Ok(result)
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    pub(crate) const ERROR: f32 = f32::EPSILON;

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
