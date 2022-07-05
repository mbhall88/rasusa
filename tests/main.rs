use assert_cmd::prelude::*;
// Add methods on commands
use predicates::prelude::*;
use std::process::Command; // Run programs // Used for writing assertions

#[test]
fn input_file_doesnt_exist() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::main_binary()?;
    cmd.args(vec!["-i", "file/doesnt/exist.fa", "-g", "5mb", "-c", "20"]);
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("does not exist"));

    Ok(())
}

#[test]
fn output_file_in_nonexistant_dir() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::main_binary()?;
    cmd.args(vec![
        "-i",
        "tests/cases/file1.fq.gz",
        "-g",
        "5mb",
        "-c",
        "20",
        "-o",
        "dir/doesnt/exists/out.fq.gz",
    ]);
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("No such file"));

    Ok(())
}

#[test]
fn valid_inputs_raises_no_errors() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::main_binary()?;
    cmd.args(vec![
        "-i",
        "tests/cases/file1.fq.gz",
        "-g",
        "5mb",
        "-c",
        "20",
    ]);

    cmd.assert().success();

    Ok(())
}

#[test]
fn input_and_output_filetypes_different_raises_no_errors() -> Result<(), Box<dyn std::error::Error>>
{
    let mut cmd = Command::main_binary()?;
    cmd.args(vec![
        "-i",
        "tests/cases/file1.fq.gz",
        "-g",
        "5mb",
        "-c",
        "20",
        "-o",
        "/tmp/out.fasta",
    ]);

    cmd.assert().success();

    Ok(())
}

#[test]
fn invalid_input_and_output_combination_raises_error() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::main_binary()?;
    cmd.args(vec![
        "-i",
        "tests/cases/file1.fq.gz",
        "-g",
        "5mb",
        "-c",
        "20",
        "-o",
        "/tmp/out.fasta",
        "-o",
        "/tmp/out2.fq",
    ]);

    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("Got 1 --input but 2 --output"));

    Ok(())
}

#[test]
fn num_instead_of_coverage_based_raises_no_errors() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::main_binary()?;
    cmd.args(vec!["-i", "tests/cases/file1.fq.gz", "-n", "5m"]);

    cmd.assert().success();

    Ok(())
}

#[test]
fn frac_instead_of_coverage_based_raises_no_errors() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::main_binary()?;
    cmd.args(vec!["-i", "tests/cases/file1.fq.gz", "-f", "0.2"]);

    cmd.assert().success();

    Ok(())
}

#[test]
fn unequal_number_of_reads_raises_error() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::main_binary()?;
    cmd.args(vec![
        "-i",
        "tests/cases/file1.fq.gz",
        "tests/cases/r2.fq.gz",
        "-g",
        "5mb",
        "-c",
        "20",
        "-o",
        "/tmp/out.fq",
        "-o",
        "/tmp/out2.fq",
    ]);

    cmd.assert().failure().stderr(predicate::str::contains(
        "Illumina files are assumed to have the same number of reads",
    ));

    Ok(())
}

#[test]
fn two_valid_illumina_inputs_suceeds() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::main_binary()?;
    cmd.args(vec![
        "-i",
        "tests/cases/r1.fq.gz",
        "tests/cases/r2.fq.gz",
        "-g",
        "4",
        "-c",
        "2",
        "-o",
        "/tmp/out.fq",
        "-o",
        "/tmp/out2.fq",
    ]);

    cmd.assert().success();

    Ok(())
}
