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

#[test]
fn reads_bam_num() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![READS, "tests/cases/ubam/single_ubam.bam", "-n", "10"]);
    cmd.assert().success();
    Ok(())
}

#[test]
fn reads_bam_frac() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![READS, "tests/cases/ubam/single_ubam.bam", "-f", "0.1"]);
    cmd.assert().success();
    Ok(())
}

#[test]
fn reads_bam_coverage() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/ubam/single_ubam.bam",
        "-c",
        "1",
        "-g",
        "10kb",
    ]);
    cmd.assert().success();
    Ok(())
}

#[test]
fn reads_bam_to_sam() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    let temp_dir = tempfile::tempdir().unwrap();
    let out = temp_dir.path().join("out.sam");
    cmd.args(vec![
        READS,
        "tests/cases/ubam/single_ubam.bam",
        "-n",
        "10",
        "-o",
        out.to_str().unwrap(),
    ]);
    cmd.assert().success();
    
    // Verify it is indeed SAM
    let content = std::fs::read_to_string(out).unwrap();
    assert!(content.starts_with("@HD"));
    Ok(())
}

#[test]
fn reads_reproducibility() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    let temp_dir1 = tempfile::tempdir().unwrap();
    let temp_dir2 = tempfile::tempdir().unwrap();
    let out1 = temp_dir1.path().join("out.bam");
    let out2 = temp_dir2.path().join("out.bam");
    let seed = "42";

    cmd.args(vec![
        READS,
        "tests/cases/ubam/single_ubam.bam",
        "-n",
        "10",
        "-s",
        seed,
        "-o",
        out1.to_str().unwrap(),
    ]);
    cmd.assert().success();

    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/ubam/single_ubam.bam",
        "-n",
        "10",
        "-s",
        seed,
        "-o",
        out2.to_str().unwrap(),
    ]);
    cmd.assert().success();

    let mut reader1 = noodles_util::alignment::io::reader::Builder::default()
        .build_from_path(out1)
        .unwrap();
    let header1 = reader1.read_header().unwrap();
    let recs1: Vec<_> = reader1.records(&header1).map(|r| r.unwrap().name().map(|n| n.to_vec())).collect();

    let mut reader2 = noodles_util::alignment::io::reader::Builder::default()
        .build_from_path(out2)
        .unwrap();
    let header2 = reader2.read_header().unwrap();
    let recs2: Vec<_> = reader2.records(&header2).map(|r| r.unwrap().name().map(|n| n.to_vec())).collect();

    assert_eq!(recs1, recs2);
    Ok(())
}

#[test]
fn reads_paired_bam() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/ubam/paired_interleave_ubam.bam",
        "-n",
        "10",
    ]);
    cmd.assert().success();
    Ok(())
}

#[test]
fn reads_single_ubam_default_output() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    let temp_dir = tempfile::tempdir().unwrap();
    let out = temp_dir.path().join("out.bam");
    cmd.args(vec![
        READS,
        "tests/cases/ubam/single_ubam.bam",
        "-n",
        "10",
        "-o",
        out.to_str().unwrap(),
    ]);
    cmd.assert().success();
    
    // Verify it is indeed BAM
    let mut reader = noodles_util::alignment::io::reader::Builder::default()
        .build_from_path(out)
        .unwrap();
    let header = reader.read_header().unwrap();
    assert!(header.programs().as_ref().get(&b"rasusa"[..]).is_some());
    Ok(())
}

#[test]
fn reads_single_ubam_fastq_output() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    let temp_dir = tempfile::tempdir().unwrap();
    let out = temp_dir.path().join("out.fq");
    cmd.args(vec![
        READS,
        "tests/cases/ubam/single_ubam.bam",
        "-n",
        "10",
        "-o",
        out.to_str().unwrap(),
    ]);
    cmd.assert().success();
    
    // Verify it is indeed FASTQ
    let content = std::fs::read_to_string(out).unwrap();
    assert!(content.starts_with("@"));
    Ok(())
}

#[test]
fn reads_mapped_bam_errors() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/test.bam",
        "-n",
        "10",
    ]);
    cmd.assert().failure().stderr(predicate::str::contains("Error: Mapped read detected, please use `rasusa aln` for aligned data"));
    Ok(())
}

#[test]
fn reads_fastq_to_bam_errors() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/r1.fq.gz",
        "-n",
        "10",
        "-o",
        "out.bam",
    ]);
    cmd.assert().failure().stderr(predicate::str::contains("Conversion from FASTA/FASTQ to Bam is not supported"));
    Ok(())
}

#[test]
fn reads_bam_to_fasta() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    let temp_dir = tempfile::tempdir().unwrap();
    let out = temp_dir.path().join("out.fasta");
    cmd.args(vec![
        READS,
        "tests/cases/ubam/single_ubam.bam",
        "-n",
        "10",
        "-o",
        out.to_str().unwrap(),
    ]);
    cmd.assert().success();

    let content = std::fs::read_to_string(out).unwrap();
    assert!(content.starts_with(">"));
    assert!(!content.contains("+"));
    Ok(())
}
