//! Command-line tests for `rasusa reads`.
//!
//! These run the built binary, so they live here as integration tests: Cargo only sets
//! `CARGO_BIN_EXE_rasusa` for integration tests, and from Rust 1.100 `assert_cmd` cannot find
//! the binary from a unit test.

use assert_cmd::Command;

const SUB: &str = "reads";

#[test]
fn no_args_given_raises_error() {
    let passed_args = vec![SUB];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().failure();
}

#[test]
fn no_input_file_given_raises_error() {
    let passed_args = vec![SUB, "-c", "30", "-g", "3mb"];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().failure();
}

#[test]
fn no_coverage_given_raises_error() {
    let infile = "tests/cases/r1.fq.gz";
    let passed_args = vec![SUB, infile, "-g", "3mb"];
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

    cmd.args(passed_args).assert().failure();
}

#[test]
fn invalid_coverage_given_raises_error() {
    let passed_args = vec![SUB, "in.fq", "-g", "3mb", "-c", "foo"];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().failure();
}

#[test]
fn no_genome_size_given_raises_error() {
    let infile = "tests/cases/r1.fq.gz";
    let passed_args = vec![SUB, infile, "-c", "5"];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().failure();
}

#[test]
fn no_genome_size_but_bases_and_coverage() {
    let infile = "tests/cases/r1.fq.gz";
    let passed_args = vec![SUB, infile, "-b", "5", "-c", "7"];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().failure();
}

#[test]
fn bases_and_coverage_and_genome_size_all_allowed() {
    let infile = "tests/cases/r1.fq.gz";
    let passed_args = vec![SUB, infile, "-b", "5", "-c", "7", "-g", "5m"];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().success();
}

#[test]
fn no_genome_size_or_coverage_given_but_bases_given() {
    let infile = "tests/cases/r1.fq.gz";
    let passed_args = vec![SUB, infile, "-b", "5"];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().success();
}

#[test]
fn no_genome_size_or_coverage_given_but_num_given() {
    let infile = "tests/cases/r1.fq.gz";
    let passed_args = vec![SUB, infile, "-n", "5"];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().success();
}

#[test]
fn no_genome_size_or_coverage_given_but_frac_given() {
    let infile = "tests/cases/r1.fq.gz";
    let passed_args = vec![SUB, infile, "-f", "5"];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().success();
}

#[test]
fn num_and_coverage_and_genome_size_not_allowed() {
    let infile = "tests/cases/r1.fq.gz";
    let passed_args = vec![SUB, infile, "-n", "5", "-c", "7", "-g", "5m"];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().failure();
}

#[test]
fn frac_and_coverage_and_genome_size_not_allowed() {
    let infile = "tests/cases/r1.fq.gz";
    let passed_args = vec![SUB, infile, "-f", "5", "-c", "7", "-g", "5m"];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().failure();
}

#[test]
fn frac_and_bases_not_allowed() {
    let infile = "tests/cases/r1.fq.gz";
    let passed_args = vec![SUB, infile, "-f", "5", "-b", "7"];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().failure();
}

#[test]
fn frac_and_num_not_allowed() {
    let infile = "tests/cases/r1.fq.gz";
    let passed_args = vec![SUB, infile, "-f", "5", "-n", "7"];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().failure();
}

#[test]
fn bases_and_num_not_allowed() {
    let infile = "tests/cases/r1.fq.gz";
    let passed_args = vec![SUB, infile, "-b", "5", "-n", "7"];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().failure();
}

#[test]
fn invalid_genome_size_given_raises_error() {
    let passed_args = vec![SUB, "in.fq", "-c", "5", "-g", "8jb"];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().failure();
}

#[test]
fn faidx_given_as_genome_size() {
    let infile = "tests/cases/r1.fq.gz";
    let faidx = "tests/cases/h37rv.fa.fai";
    let passed_args = vec![SUB, infile, "-c", "5", "-g", faidx];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().success();
}

#[test]
fn invalid_seed_given_raises_error() {
    let infile = "tests/cases/r1.fq.gz";
    let passed_args = vec![SUB, infile, "-c", "5", "-g", "8mb", "-s", "foo"];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().failure();
}

#[test]
fn all_valid_args_parsed_as_expected() {
    let infile = "tests/cases/r1.fq.gz";
    let passed_args = vec![
        SUB,
        infile,
        "-c",
        "5",
        "-g",
        "8mb",
        "-s",
        "88",
        "-o",
        "/dev/null",
    ];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().success();
}

#[test]
fn two_outputs_passed_with_one_option_flag() {
    let infile = "tests/cases/r1.fq.gz";
    let passed_args = vec![
        SUB,
        infile,
        infile,
        "-c",
        "5",
        "-g",
        "8mb",
        "-s",
        "88",
        "-o",
        "/dev/null",
        "/dev/null",
    ];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().failure();
}

#[test]
fn three_inputs_raises_error() {
    let infile = "tests/cases/r1.fq.gz";
    let passed_args = vec![
        SUB,
        infile,
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

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().failure();
}

#[test]
fn three_outputs_raises_error() {
    let infile = "tests/cases/r1.fq.gz";
    let passed_args = vec![
        SUB,
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
        "-o",
        "out.fq",
        "-o",
        "out.fq",
    ];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().failure();
}

#[test]
fn one_input_no_output_is_ok() {
    let infile = "tests/cases/r1.fq.gz";
    let passed_args = vec![SUB, infile, "-c", "5", "-g", "8mb", "-s", "88"];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().success();
}

#[test]
fn two_inputs_one_output_raises_error() {
    let infile = "tests/cases/r1.fq.gz";
    let passed_args = vec![
        SUB, infile, infile, "-c", "5", "-g", "8mb", "-s", "88", "-o", "out.fq",
    ];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().failure();
}

#[test]
fn one_input_two_outputs_raises_error() {
    let infile = "tests/cases/r1.fq.gz";
    let passed_args = vec![
        SUB, infile, "-c", "5", "-g", "8mb", "-s", "88", "-o", "out.fq", "-o", "out.fq",
    ];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().failure();
}

#[test]
fn two_input_two_outputs_is_ok() {
    let infile = "tests/cases/r1.fq.gz";
    let passed_args = vec![
        SUB, infile, infile, "-c", "5", "-g", "8mb", "-s", "88", "-o", "out.fq", "-o", "out.fq",
    ];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().success();
}

#[test]
fn two_input_two_outputs_is_ok_when_positional_args_at_end() {
    let infile = "tests/cases/r1.fq.gz";
    let passed_args = vec![
        SUB, "-c", "5", "-g", "8mb", "-s", "88", "-o", "out.fq", "-o", "out.fq", infile, infile,
    ];

    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();
    cmd.args(passed_args).assert().success();
}
