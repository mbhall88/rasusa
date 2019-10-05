use regex::Regex;
use snafu::Snafu;
use std::ops::Mul;
use std::path::PathBuf;
use std::str::FromStr;
use structopt::StructOpt;

/// Randomly sub-sample reads to a specified coverage
#[derive(Debug, StructOpt)]
#[structopt()]
pub struct Cli {
    /// The fast{a,q} file to sub-sample.
    #[structopt(short = "i", long = "input", parse(from_os_str))]
    input: PathBuf,

    /// Output file, stdout if not present.
    #[structopt(short = "o", long = "output", parse(from_os_str))]
    output: Option<PathBuf>,

    /// Size of the genome to calculate coverage with respect to. i.e 4.3kb, 7Tb, 9000, 4.1MB etc.
    #[structopt(short = "g", long = "genome-size")]
    genome_size: GenomeSize,

    /// The desired coverage to sub-sample the reads to.
    #[structopt(short = "c", long = "coverage")]
    coverage: Coverage,

    /// Random seed to use.
    #[structopt(short = "s", long = "seed")]
    seed: Option<u64>,

    /// Switch on verbosity
    #[structopt(short)]
    verbose: bool,
}

#[derive(Debug, Snafu, PartialEq)]
enum Invalid {
    #[snafu(display("{} is not a valid metric suffix", suffix))]
    MetricSuffixString { suffix: String },
    #[snafu(display(
        "{} is not a valid genome size. Valid forms include 4gb, 3000, 8.7Kb etc.",
        genome_size
    ))]
    GenomeSizeString { genome_size: String },
    #[snafu(display("{} is not a valid coverage string. Coverage must be either an integer or a float and can end with an optional 'x' character", coverage))]
    CoverageValue { coverage: String },
}

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

#[derive(Debug, PartialOrd, PartialEq)]
struct GenomeSize(u64);

impl PartialEq<u64> for GenomeSize {
    fn eq(&self, other: &u64) -> bool {
        self.0 == *other
    }
}

impl FromStr for GenomeSize {
    type Err = Invalid;

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
        let metric_suffix = match captures.name("sfx") {
            Some(suffix) => MetricSuffix::from_str(suffix.as_str())?,
            None => {
                return Err(Invalid::MetricSuffixString {
                    suffix: "".to_string(),
                });
            }
        };
        Ok(GenomeSize((size * metric_suffix) as u64))
    }
}

#[derive(Debug, PartialOrd, PartialEq)]
struct Coverage(u32);

impl PartialEq<u32> for Coverage {
    fn eq(&self, other: &u32) -> bool {
        self.0 == *other
    }
}

impl FromStr for Coverage {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let re = Regex::new(r"^(?P<covg>[0-9]*\.?[0-9]+)(?i)x?$").unwrap();
        let captures = match re.captures(s) {
            Some(cap) => cap,
            None => {
                return Err(
                    "Coverage should be an int or float and can end in an x character.".to_string(),
                );
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
        let actual = GenomeSize::from_str(".77uB");

        assert!(actual.is_err());
    }

    #[test]
    fn empty_string_returns_error() {
        let actual = GenomeSize::from_str("");

        assert!(actual.is_err());
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
        let actual = Coverage::from_str("");

        assert!(actual.is_err())
    }

    #[test]
    fn non_number_coverage_returns_err() {
        let actual = Coverage::from_str("foo");

        assert!(actual.is_err())
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
}
