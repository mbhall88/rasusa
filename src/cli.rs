extern crate clap;

use clap::{App, Arg, ArgMatches};
use regex::Regex;
use std::ffi::OsString;
use std::ops::Mul;
use std::str::FromStr;

#[derive(PartialEq, Debug)]
enum MetricSuffix {
    Base,
    Kilo,
    Mega,
    Giga,
    Tera,
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

#[derive(Debug)]
struct GenomeSize(u64);

impl PartialEq<u64> for GenomeSize {
    fn eq(&self, other: &u64) -> bool {
        self.0 == *other
    }
}

impl FromStr for GenomeSize {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let text = s.to_lowercase();
        let re = Regex::new(r"(?P<size>[0-9]*\.?[0-9]+)(?P<sfx>\w*)$").unwrap();
        let captures = match re.captures(text.as_str()) {
            Some(cap) => cap,
            None => return Err(
                "Invalid genome size format. Accepted formats include: 4000, 4mb, 7Gb, 6.7Kb etc."
                    .to_string(),
            ),
        };
        let size = captures
            .name("size")
            .unwrap()
            .as_str()
            .parse::<f64>()
            .unwrap();
        let metric_suffix = match captures.name("sfx") {
            Some(x) if "b".contains(x.as_str()) => MetricSuffix::Base,
            Some(x) if "kb".contains(x.as_str()) => MetricSuffix::Kilo,
            Some(x) if "mb".contains(x.as_str()) => MetricSuffix::Mega,
            Some(x) if "gb".contains(x.as_str()) => MetricSuffix::Giga,
            Some(x) if "tb".contains(x.as_str()) => MetricSuffix::Tera,
            None => MetricSuffix::Base,
            _ => return Err("Invalid metric suffix.".to_string()),
        };
        Ok(GenomeSize((size * metric_suffix) as u64))
    }
}

#[derive(Debug)]
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

pub fn parse_args<'a, I, T>(args: I) -> ArgMatches<'a>
where
    I: IntoIterator<Item = T>,
    T: Into<OsString> + Clone,
{
    App::new("rasusa")
        .version("0.1.0")
        .author("Michael B. Hall <michael@mbh.sh>")
        .about("Randomly sub-sample reads to a specified coverage")
        .arg(
            Arg::with_name("input")
                .short("i")
                .long("input")
                .value_name("INPUT")
                .help("The fastq file to sub-sample.")
                .required(true)
                .takes_value(true)
                .display_order(1),
        )
        .arg(
            Arg::with_name("genome_size")
                .short("g")
                .long("genome-size")
                .value_name("GENOME SIZE")
                .help(
                    // todo: fix formatting of text in terminal
                    "Size of the genome to calculate coverage with respect to. \
                     This can be either an integer or an integer/float with a metric suffix. \
                     An Example of this is 4.3kb will indicate a genome size of 4300. \
                     Other examples include 7Tb, 9.1gb, 90MB, 6009 etc.",
                )
                .required(true)
                .takes_value(true),
        )
        .arg(
            Arg::with_name("coverage")
                .short("c")
                .long("coverage")
                .value_name("COVG")
                .help("The desired coverage to sub-sample the reads to.")
                .required(true)
                .takes_value(true),
        )
        .arg(
            Arg::with_name("v")
                .short("v")
                .multiple(true)
                .help("Sets the level of verbosity"),
        )
        .get_matches_from(args)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn float_multiply_with_base_unit() {
        let actual = 4.5 * MetricSuffix::Base;
        let expected = 4.5;

        assert_eq!(actual, expected)
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
}
