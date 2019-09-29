extern crate clap;

use clap::{App, Arg, ArgMatches};
use regex::Regex;
use std::ffi::OsString;

#[derive(PartialEq, Debug)]
enum MetricSuffix {
    Base,
    Kilo,
    Mega,
    Giga,
    Tera,
}

fn multiply_size_by_metric_suffix(
    size: f64,
    metric_suffix: Result<MetricSuffix, String>,
) -> Result<f64, String> {
    match metric_suffix {
        Ok(MetricSuffix::Base) => Ok(size),
        Ok(MetricSuffix::Kilo) => Ok(size * 1_000.0),
        Ok(MetricSuffix::Mega) => Ok(size * 1_000_000.0),
        Ok(MetricSuffix::Giga) => Ok(size * 1_000_000_000.0),
        Ok(MetricSuffix::Tera) => Ok(size * 1_000_000_000_000.0),
        Err(msg) => Err(msg),
    }
}

fn parse_genome_size(s: String) -> Result<u64, String> {
    let text = s.to_lowercase();
    let re = Regex::new(r"(?P<size>[0-9]*\.?[0-9]+)(?P<sfx>\w*)$").unwrap();
    let captures = match re.captures(text.as_str()) {
        Some(cap) => cap,
        None => return Err("No matches found in genome size string".to_string()),
    };
    let size = match captures.name("size") {
        Some(x) => x.as_str().parse::<f64>().unwrap(),
        None => 0 as f64,
    };
    let metric_suffix = match captures.name("sfx") {
        Some(x) if "b".contains(x.as_str()) => Ok(MetricSuffix::Base),
        Some(x) if "kb".contains(x.as_str()) => Ok(MetricSuffix::Kilo),
        Some(x) if "mb".contains(x.as_str()) => Ok(MetricSuffix::Mega),
        Some(x) if "gb".contains(x.as_str()) => Ok(MetricSuffix::Giga),
        Some(x) if "tb".contains(x.as_str()) => Ok(MetricSuffix::Tera),
        None => Ok(MetricSuffix::Base),
        _ => Err("Invalid metric suffix.".to_string()),
    };
    match multiply_size_by_metric_suffix(size, metric_suffix) {
        Ok(x) => Ok(x as u64),
        Err(msg) => Err(msg),
    }
}

fn is_genome_size_valid(arg: String) -> Result<(), String> {
    match parse_genome_size(arg) {
        Ok(_) => Ok(()),
        Err(msg) => Err(msg),
    }
}

fn parse_coverage(arg: String) -> Result<u32, String> {
    let re = Regex::new(r"^(?P<covg>[0-9]*\.?[0-9]+)(?i)x?$").unwrap();
    let captures = match re.captures(arg.as_str()) {
        Some(cap) => cap,
        None => {
            return Err(
                "Coverage should be an int or float and can end in an x character.".to_string(),
            );
        }
    };
    match captures.name("covg") {
        Some(x) => Ok(x.as_str().parse::<f64>().unwrap() as u32),
        None => {
            Err("Coverage should be an int or float and can end in an x character.".to_string())
        }
    }
}

fn is_coverage_valid(arg: String) -> Result<(), String> {
    match parse_coverage(arg) {
        Ok(0) => Err("Coverage of 0 given.".to_string()),
        Ok(_) => Ok(()),
        Err(msg) => Err(msg),
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
                .takes_value(true)
                .validator(is_genome_size_valid),
        )
        .arg(
            Arg::with_name("coverage")
                .short("c")
                .long("coverage")
                .value_name("COVG")
                .help("The desired coverage to sub-sample the reads to.")
                .required(true)
                .takes_value(true)
                .validator(is_coverage_valid),
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
    fn integer_only_returns_integer() {
        let arg = String::from("6");

        let actual = parse_genome_size(arg).unwrap();
        let expected = 6;

        assert_eq!(actual, expected);
    }

    #[test]
    fn float_only_returns_integer() {
        let arg = String::from("6.5");

        let actual = parse_genome_size(arg).unwrap();
        let expected = 6;

        assert_eq!(actual, expected);
    }

    #[test]
    fn int_and_suffix_returns_multiplied_int() {
        let arg = String::from("5mb");

        let actual = parse_genome_size(arg).unwrap();
        let expected = 5_000_000 as u64;

        assert_eq!(actual, expected);
    }

    #[test]
    fn float_and_suffix_returns_multiplied_float_as_int() {
        let arg = String::from("5.4kB");

        let actual = parse_genome_size(arg).unwrap();
        let expected = 5_400 as u64;

        assert_eq!(actual, expected);
    }

    #[test]
    fn float_without_leading_int_and_suffix_returns_multiplied_float_as_int() {
        let arg = String::from(".77G");

        let actual = parse_genome_size(arg).unwrap();
        let expected = 770_000_000 as u64;

        assert_eq!(actual, expected);
    }

    #[test]
    fn int_and_tera_suffix_returns_multiplied_int() {
        let arg = String::from("7TB");

        let actual = parse_genome_size(arg).unwrap();
        let expected = 7_000_000_000_000 as u64;

        assert_eq!(actual, expected);
    }

    #[test]
    fn int_and_base_suffix_returns_int_without_scaling() {
        let arg = String::from("7B");

        let actual = parse_genome_size(arg).unwrap();
        let expected = 7 as u64;

        assert_eq!(actual, expected);
    }

    #[test]
    fn invalid_suffix_returns_err() {
        let arg = String::from(".77uB");

        let actual = parse_genome_size(arg);

        assert!(actual.is_err());
    }

    #[test]
    fn empty_string_returns_error() {
        let arg = String::from("");

        let actual = parse_genome_size(arg);

        assert!(actual.is_err());
    }

    #[test]
    fn suffix_with_no_size_returns_error() {
        let arg = String::from("gb");

        let actual = parse_genome_size(arg);

        assert!(actual.is_err());
    }

    #[test]
    fn empty_string_is_invalid_genome_size() {
        let arg = String::from("");

        assert!(is_genome_size_valid(arg).is_err())
    }

    #[test]
    fn int_only_is_valid() {
        let arg = String::from("789");

        assert!(is_genome_size_valid(arg).is_ok())
    }

    #[test]
    fn float_only_is_valid() {
        let arg = String::from("79.6");

        assert!(is_genome_size_valid(arg).is_ok())
    }

    #[test]
    fn int_with_valid_suffix_is_valid() {
        let arg = String::from("79g");

        assert!(is_genome_size_valid(arg).is_ok())
    }

    #[test]
    fn int_with_invalid_suffix_is_invalid() {
        let arg = String::from("79p");

        assert!(is_genome_size_valid(arg).is_err())
    }

    #[test]
    fn float_with_valid_suffix_is_valid() {
        let arg = String::from("79.765t");

        assert!(is_genome_size_valid(arg).is_ok())
    }

    #[test]
    fn float_with_invalid_suffix_is_invalid() {
        let arg = String::from("79.765eb");

        assert!(is_genome_size_valid(arg).is_err())
    }

    #[test]
    fn int_coverage_returns_int() {
        let arg = String::from("56");

        let actual = parse_coverage(arg).unwrap();
        let expected = 56 as u32;

        assert_eq!(actual, expected)
    }

    #[test]
    fn float_coverage_returns_int() {
        let arg = String::from("56.6");

        let actual = parse_coverage(arg).unwrap();
        let expected = 56 as u32;

        assert_eq!(actual, expected)
    }

    #[test]
    fn empty_coverage_returns_err() {
        let arg = String::from("");

        let actual = parse_coverage(arg);

        assert!(actual.is_err())
    }

    #[test]
    fn non_number_coverage_returns_err() {
        let arg = String::from("foo");

        let actual = parse_coverage(arg);

        assert!(actual.is_err())
    }

    #[test]
    fn zero_coverage_returns_zero() {
        let arg = String::from("0");

        let actual = parse_coverage(arg).unwrap();
        let expected = 0 as u32;

        assert_eq!(actual, expected)
    }

    #[test]
    fn int_ending_in_x_coverage_returns_int() {
        let arg = String::from("1X");

        let actual = parse_coverage(arg).unwrap();
        let expected = 1 as u32;

        assert_eq!(actual, expected)
    }

    #[test]
    fn float_ending_in_x_coverage_returns_int() {
        let arg = String::from("1.9x");

        let actual = parse_coverage(arg).unwrap();
        let expected = 1 as u32;

        assert_eq!(actual, expected)
    }

    #[test]
    fn empty_string_is_invalid_coverage_string() {
        let covg = String::from("");

        assert!(is_coverage_valid(covg).is_err())
    }

    #[test]
    fn zero_is_invalid_coverage_string() {
        let covg = String::from("0");

        assert!(is_coverage_valid(covg).is_err())
    }

    #[test]
    fn int_is_valid_coverage_string() {
        let covg = String::from("90");

        assert!(is_coverage_valid(covg).is_ok())
    }

    #[test]
    fn float_is_valid_coverage_string() {
        let covg = String::from("90.9");

        assert!(is_coverage_valid(covg).is_ok())
    }

    #[test]
    fn word_is_invalid_coverage_string() {
        let covg = String::from("foo");

        assert!(is_coverage_valid(covg).is_err())
    }

    #[test]
    fn int_ending_in_x_is_valid_coverage_string() {
        let covg = String::from("50x");

        assert!(is_coverage_valid(covg).is_ok())
    }

    #[test]
    fn float_ending_in_x_is_valid_coverage_string() {
        let covg = String::from("50.8X");

        assert!(is_coverage_valid(covg).is_ok())
    }
}
