use assert_cmd::prelude::*; // Add methods on commands
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
