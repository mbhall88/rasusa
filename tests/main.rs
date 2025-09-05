use assert_cmd::Command;
// Add methods on commands
use predicates::prelude::*;

const BIN: &str = "rasusa";
const READS: &str = "reads";
#[test]
fn input_file_doesnt_exist() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![READS, "file/doesnt/exist.fa", "-g", "5mb", "-c", "20"]);
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("does not exist"));

    Ok(())
}

#[test]
fn output_file_in_nonexistant_dir() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
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
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
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
fn valid_inputs_but_strict_raises_error_coverage() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/file1.fq.gz",
        "-g",
        "5mb",
        "-c",
        "20",
        "--strict",
    ]);

    cmd.assert().failure().stderr(predicate::str::contains(
        "is not possible as the actual coverage",
    ));

    Ok(())
}

#[test]
fn valid_inputs_but_strict_raises_error_bases() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/file1.fq.gz",
        "-b",
        "5mb",
        "--strict",
    ]);

    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("is more than the input"));

    Ok(())
}

#[test]
fn valid_inputs_but_strict_raises_error_num_reads() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/file1.fq.gz",
        "-n",
        "500",
        "--strict",
    ]);

    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("is more than the input"));

    Ok(())
}

#[test]
fn valid_inputs_but_strict_raises_error_frac() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/file1.fq.gz",
        "-f",
        "0.05",
        "--strict",
    ]);

    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("was rounded to 0"));

    Ok(())
}

#[test]
fn input_and_output_filetypes_different_raises_no_errors() -> Result<(), Box<dyn std::error::Error>>
{
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
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
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
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
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![READS, "tests/cases/file1.fq.gz", "-n", "5m"]);

    cmd.assert().success();

    Ok(())
}

#[test]
fn frac_instead_of_coverage_based_raises_no_errors() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![READS, "tests/cases/file1.fq.gz", "-f", "0.2"]);

    cmd.assert().success();

    Ok(())
}

#[test]
fn unequal_number_of_reads_raises_error() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
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
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
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

#[test]
fn num_from_each_with_paired_reads() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/r1.fq.gz",
        "tests/cases/r2.fq.gz",
        "-n",
        "1",
        "-o",
        "/tmp/out.fq",
        "-o",
        "/tmp/out2.fq",
    ]);

    cmd.assert().success();

    Ok(())
}

#[test]
fn frac_from_each_with_paired_reads() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/r1.fq.gz",
        "tests/cases/r2.fq.gz",
        "-f",
        "10",
        "-o",
        "/tmp/out.fq",
        "-o",
        "/tmp/out2.fq",
    ]);

    cmd.assert().success();

    Ok(())
}
