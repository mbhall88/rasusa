//! Command-line tests for `rasusa aln`.
//!
//! These run the built binary, so they live here as integration tests: Cargo only sets
//! `CARGO_BIN_EXE_rasusa` for integration tests, and from Rust 1.100 `assert_cmd` cannot find
//! the binary from a unit test.

use assert_cmd::Command;

const SUB: &str = "aln";

#[test]
fn no_coverage_given_raises_error() {
    let infile = "tests/cases/test.bam";
    let passed_args = vec![SUB, infile];
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

    cmd.args(passed_args).assert().failure();
}

#[test]
fn zero_coverage_raises_error_stream() {
    let infile = "tests/cases/test.bam";
    let passed_args = vec![SUB, infile, "-c", "0"];
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

    cmd.args(passed_args).assert().failure();
}

#[test]
fn zero_coverage_raises_error_fetch() {
    let infile = "tests/cases/test.bam";
    let passed_args = vec![
        SUB,
        infile,
        "-c",
        "0",
        "--strategy",
        "fetch",
        "--step-size",
        "5000",
    ];
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

    cmd.args(passed_args).assert().failure();
}

#[test]
fn bam_with_regions_of_zero_coverage_doesnt_endless_loop_stream() {
    let infile = "tests/cases/test.bam";
    let passed_args = vec![SUB, infile, "-c", "1"];
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

    cmd.args(passed_args).assert().success();
}

#[test]
fn bam_with_regions_of_zero_coverage_doesnt_endless_loop_fetch() {
    let infile = "tests/cases/test.bam";
    let passed_args = vec![SUB, infile, "-c", "1", "--strategy", "fetch"];
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

    cmd.args(passed_args).assert().success();
}

#[test]
fn excess_coverage_doesnt_endless_loop_stream() {
    let infile = "tests/cases/test.bam";
    let passed_args = vec![SUB, infile, "-c", "10000"];
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

    cmd.args(passed_args).assert().success();
}

#[test]
fn excess_coverage_doesnt_endless_loop_fetch() {
    let infile = "tests/cases/test.bam";
    let passed_args = vec![SUB, infile, "-c", "10000", "--strategy", "fetch"];
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

    cmd.args(passed_args).assert().success();
}

#[test]
fn bam_no_index_is_ok_stream() {
    let infile = "tests/cases/no_index.bam";
    let passed_args = vec![SUB, infile, "-c", "1"];
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

    cmd.args(passed_args).assert().success();
}

#[test]
fn bam_no_index_fails_fetch() {
    let infile = "tests/cases/no_index.bam";
    let passed_args = vec![SUB, infile, "-c", "1", "--strategy", "fetch"];
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

    cmd.args(passed_args).assert().failure();
}

#[test]
fn bam_with_no_start_or_end_regions_and_missing_chromosomes() {
    let infile = "tests/cases/no_start_end.bam";
    let passed_args = vec![SUB, infile, "-c", "1"];
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

    cmd.args(passed_args).assert().success();
}

#[test]
fn bam_is_not_sorted_fails() {
    let infile = "tests/cases/test_not_sorted.bam";
    let passed_args = vec![SUB, infile, "-c", "1"];
    let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

    cmd.args(passed_args).assert().failure();
}
