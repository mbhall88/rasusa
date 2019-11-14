use regex::Regex;
use snafu::Snafu;
use std::ops::{Div, Mul};
use std::path::PathBuf;
use std::str::FromStr;
use structopt::StructOpt;

/// Randomly subsample reads to a specified coverage.
#[derive(Debug, StructOpt)]
#[structopt()]
pub struct Cli {
    /// The fast{a,q} file to sub-sample.
    #[structopt(short = "i", long = "input", parse(from_os_str))]
    pub input: PathBuf,

    /// Output file, stdout if not present.
    #[structopt(short = "o", long = "output", parse(from_os_str))]
    pub output: Option<PathBuf>,

    /// Size of the genome to calculate coverage with respect to. i.e 4.3kb, 7Tb, 9000, 4.1MB etc.
    #[structopt(short = "g", long = "genome-size")]
    pub genome_size: GenomeSize,

    /// The desired coverage to sub-sample the reads to.
    #[structopt(short = "c", long = "coverage")]
    pub coverage: Coverage,

    /// Random seed to use.
    #[structopt(short = "s", long = "seed")]
    pub seed: Option<u64>,

    /// Switch on verbosity.
    #[structopt(short)]
    pub verbose: bool,
}

/// A collection of custom errors relating to the command line interface for this package.
#[derive(Debug, Snafu, PartialEq)]
pub enum Invalid {
    /// Indicates that a string cannot be parsed into a [`MetricSuffix`](#metricsuffix).
    #[snafu(display("{} is not a valid metric suffix", suffix))]
    MetricSuffixString { suffix: String },

    /// Indicates that a string cannot be parsed into a [`GenomeSize`](#genomesize).
    #[snafu(display(
        "{} is not a valid genome size. Valid forms include 4gb, 3000, 8.7Kb etc.",
        genome_size
    ))]
    GenomeSizeString { genome_size: String },

    /// Indicates that a string cannot be parsed into a [`Coverage`](#coverage).
    #[snafu(display("{} is not a valid coverage string. Coverage must be either an integer or a float and can end with an optional 'x' character", coverage))]
    CoverageValue { coverage: String },
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
    type Err = Invalid;

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
                return Err(Invalid::MetricSuffixString {
                    suffix: suffix.to_string(),
                });
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

impl FromStr for GenomeSize {
    type Err = Invalid;

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
            None => {
                return Err(Invalid::GenomeSizeString {
                    genome_size: s.to_string(),
                });
            }
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
        self.0 * u64::from(rhs.0)
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
pub struct Coverage(u32);

/// Allow for comparison of a `u32` and a `Coverage`.
///
/// # Example
///
/// ```rust
/// assert!(Coverage(10) == 10)
/// ```
impl PartialEq<u32> for Coverage {
    fn eq(&self, other: &u32) -> bool {
        self.0 == *other
    }
}

impl FromStr for Coverage {
    type Err = Invalid;

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
                return Err(Invalid::CoverageValue {
                    coverage: s.to_string(),
                });
            }
        };
        Ok(Coverage(
            captures
                .name("covg")
                .unwrap()
                .as_str()
                .parse::<f64>()
                .unwrap() as u32,
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
        u64::from(self.0) * rhs.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_args_given_raises_error() {
        let passed_args = vec!["rasusa"];
        let args: Result<Cli, clap::Error> = Cli::from_iter_safe(passed_args);

        let actual = args.unwrap_err().kind;
        let expected = clap::ErrorKind::MissingRequiredArgument;

        assert_eq!(actual, expected)
    }

    #[test]
    fn no_input_file_given_raises_error() {
        let passed_args = vec!["rasusa", "-c", "30", "-g", "3mb"];
        let args: Result<Cli, clap::Error> = Cli::from_iter_safe(passed_args);

        let actual = args.unwrap_err().kind;
        let expected = clap::ErrorKind::MissingRequiredArgument;

        assert_eq!(actual, expected)
    }

    #[test]
    fn no_coverage_given_raises_error() {
        let passed_args = vec!["rasusa", "-i", "in.fq", "-g", "3mb"];
        let args: Result<Cli, clap::Error> = Cli::from_iter_safe(passed_args);

        let actual = args.unwrap_err().kind;
        let expected = clap::ErrorKind::MissingRequiredArgument;

        assert_eq!(actual, expected)
    }

    #[test]
    fn invalid_coverage_given_raises_error() {
        let passed_args = vec!["rasusa", "-i", "in.fq", "-g", "3mb", "-c", "foo"];
        let args: Result<Cli, clap::Error> = Cli::from_iter_safe(passed_args);

        let actual = args.unwrap_err().kind;
        let expected = clap::ErrorKind::ValueValidation;

        assert_eq!(actual, expected)
    }

    #[test]
    fn no_genome_size_given_raises_error() {
        let passed_args = vec!["rasusa", "-i", "in.fq", "-c", "5"];
        let args: Result<Cli, clap::Error> = Cli::from_iter_safe(passed_args);

        let actual = args.unwrap_err().kind;
        let expected = clap::ErrorKind::MissingRequiredArgument;

        assert_eq!(actual, expected)
    }

    #[test]
    fn invalid_genome_size_given_raises_error() {
        let passed_args = vec!["rasusa", "-i", "in.fq", "-c", "5", "-g", "8jb"];
        let args: Result<Cli, clap::Error> = Cli::from_iter_safe(passed_args);

        let actual = args.unwrap_err().kind;
        let expected = clap::ErrorKind::ValueValidation;

        assert_eq!(actual, expected)
    }

    #[test]
    fn invalid_seed_given_raises_error() {
        let passed_args = vec!["rasusa", "-i", "in.fq", "-c", "5", "-g", "8mb", "-s", "foo"];
        let args: Result<Cli, clap::Error> = Cli::from_iter_safe(passed_args);

        let actual = args.unwrap_err().kind;
        let expected = clap::ErrorKind::ValueValidation;

        assert_eq!(actual, expected)
    }

    #[test]
    fn all_valid_args_parsed_as_expected() {
        let passed_args = vec![
            "rasusa",
            "-i",
            "in.fq",
            "-c",
            "5",
            "-g",
            "8mb",
            "-s",
            "88",
            "-o",
            "my/output/file.fq",
        ];
        let args = Cli::from_iter_safe(passed_args).unwrap();

        assert_eq!(args.input, PathBuf::from_str("in.fq").unwrap());
        assert_eq!(args.coverage, Coverage(5));
        assert_eq!(args.genome_size, GenomeSize(8_000_000));
        assert_eq!(args.seed, Some(88));
        assert_eq!(args.output, PathBuf::from_str("my/output/file.fq").ok())
    }

    #[test]
    fn float_multiply_with_base_unit() {
        let actual = 4.5 * MetricSuffix::Base;
        let expected = 4.5;

        let diff = (actual.abs() - expected).abs();
        assert!(diff < std::f64::EPSILON)
    }

    #[test]
    fn integer_only_returns_integer() {
        let actual = GenomeSize::from_str("6").unwrap();
        let expected = 6 as u64;

        assert_eq!(actual, expected);
    }

    #[test]
    fn float_only_returns_integer() {
        let actual = GenomeSize::from_str("6.5").unwrap();
        let expected = 6 as u64;

        assert_eq!(actual, expected);
    }

    #[test]
    fn int_and_suffix_returns_multiplied_int() {
        let actual = GenomeSize::from_str("5mb").unwrap();
        let expected = 5_000_000 as u64;

        assert_eq!(actual, expected);
    }

    #[test]
    fn float_and_suffix_returns_multiplied_float_as_int() {
        let actual = GenomeSize::from_str("5.4kB").unwrap();
        let expected = 5_400 as u64;

        assert_eq!(actual, expected);
    }

    #[test]
    fn float_without_leading_int_and_suffix_returns_multiplied_float_as_int() {
        let actual = GenomeSize::from_str(".77G").unwrap();
        let expected = 770_000_000 as u64;

        assert_eq!(actual, expected);
    }

    #[test]
    fn int_and_tera_suffix_returns_multiplied_int() {
        let actual = GenomeSize::from_str("7TB").unwrap();
        let expected = 7_000_000_000_000 as u64;

        assert_eq!(actual, expected);
    }

    #[test]
    fn int_and_base_suffix_returns_int_without_scaling() {
        let actual = GenomeSize::from_str("7B").unwrap();
        let expected = 7 as u64;

        assert_eq!(actual, expected);
    }

    #[test]
    fn invalid_suffix_returns_err() {
        let genome_size = String::from(".77uB");
        let actual = GenomeSize::from_str(genome_size.as_str()).unwrap_err();
        let expected = Invalid::MetricSuffixString {
            suffix: String::from("ub"),
        };

        assert_eq!(actual, expected);
    }

    #[test]
    fn empty_string_returns_error() {
        let actual = GenomeSize::from_str("").unwrap_err();
        let expected = Invalid::GenomeSizeString {
            genome_size: String::from(""),
        };

        assert_eq!(actual, expected);
    }

    #[test]
    fn suffix_with_no_size_returns_error() {
        let actual = GenomeSize::from_str("gb");

        assert!(actual.is_err());
    }

    #[test]
    fn int_coverage_returns_int() {
        let actual = Coverage::from_str("56").unwrap();
        let expected = 56 as u32;

        assert_eq!(actual, expected)
    }

    #[test]
    fn float_coverage_returns_int() {
        let actual = Coverage::from_str("56.6").unwrap();
        let expected = 56 as u32;

        assert_eq!(actual, expected)
    }

    #[test]
    fn empty_coverage_returns_err() {
        let coverage = String::from("");

        let actual = Coverage::from_str(coverage.as_str()).unwrap_err();
        let expected = Invalid::CoverageValue { coverage };

        assert_eq!(actual, expected)
    }

    #[test]
    fn non_number_coverage_returns_err() {
        let coverage = String::from("foo");

        let actual = Coverage::from_str(coverage.as_str()).unwrap_err();
        let expected = Invalid::CoverageValue { coverage };

        assert_eq!(actual, expected)
    }

    #[test]
    fn zero_coverage_returns_zero() {
        let actual = Coverage::from_str("0").unwrap();
        let expected = 0 as u32;

        assert_eq!(actual, expected)
    }

    #[test]
    fn int_ending_in_x_coverage_returns_int() {
        let actual = Coverage::from_str("1X").unwrap();
        let expected = 1 as u32;

        assert_eq!(actual, expected)
    }

    #[test]
    fn float_ending_in_x_coverage_returns_int() {
        let actual = Coverage::from_str("1.9X").unwrap();
        let expected = 1 as u32;

        assert_eq!(actual, expected)
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
        let expected = Invalid::MetricSuffixString { suffix };

        assert_eq!(actual, expected)
    }

    #[test]
    fn multiply_genome_size_by_coverage() {
        let genome_size = GenomeSize::from_str("4.2kb").unwrap();
        let covg = Coverage::from_str("10").unwrap();

        let actual = genome_size * covg;
        let expected: u64 = 42_000;

        assert_eq!(actual, expected)
    }

    #[test]
    fn multiply_coverage_by_genome_size() {
        let genome_size = GenomeSize::from_str("4.2kb").unwrap();
        let covg = Coverage::from_str("10").unwrap();

        let actual = covg * genome_size;
        let expected: u64 = 42_000;

        assert_eq!(actual, expected)
    }

    #[test]
    fn divide_u64_by_genome_size() {
        let x: u64 = 210;
        let size = GenomeSize(200);

        let actual = x / size;
        let expected = 1.05;

        let diff = (actual.abs() - expected).abs();
        assert!(diff < std::f64::EPSILON)
    }
}
