use assert_cmd::prelude::*;
// Add methods on commands
use predicates::prelude::*;
use std::process::Command; // Run programs // Used for writing assertions

#[test]
fn invalid_input_file() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::main_binary()?;
    cmd.args(vec!["-i file/doesnt/exist.fx", "-g", "5mb", "-c", "20"]);
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("is not fasta or fastq"));

    Ok(())
}

#[test]
fn input_file_doesnt_exist() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::main_binary()?;
    cmd.args(vec!["-i", "file/doesnt/exist.fa", "-g", "5mb", "-c", "20"]);
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("No such file"));

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
fn input_and_output_filetypes_different_raises_no_errors_but_warns(
) -> Result<(), Box<dyn std::error::Error>> {
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

    cmd.assert()
        .success()
        .stderr(predicate::str::contains("file types are not the same"));

    Ok(())
}
