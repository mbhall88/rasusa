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
    let recs1: Vec<_> = reader1
        .records(&header1)
        .map(|r| r.unwrap().name().map(|n| n.to_vec()))
        .collect();

    let mut reader2 = noodles_util::alignment::io::reader::Builder::default()
        .build_from_path(out2)
        .unwrap();
    let header2 = reader2.read_header().unwrap();
    let recs2: Vec<_> = reader2
        .records(&header2)
        .map(|r| r.unwrap().name().map(|n| n.to_vec()))
        .collect();

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
    cmd.args(vec![READS, "tests/cases/test.bam", "-n", "10"]);
    cmd.assert().failure().stderr(predicate::str::contains(
        "Error: Mapped read detected, please use `rasusa aln` for aligned data",
    ));
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
    cmd.assert().failure().stderr(predicate::str::contains(
        "Conversion from FASTA/FASTQ to Bam is not supported",
    ));
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
use std::fs;

#[test]
fn reads_format_combinations() -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempfile::tempdir().unwrap();
    let inputs = vec![
        "tests/cases/ubam/test.fasta",
        "tests/cases/ubam/test.fastq",
        "tests/cases/ubam/single_usam.sam",
        "tests/cases/ubam/single_ubam.bam",
        "tests/cases/ubam/single_ucram.cram",
    ];

    let outputs = vec!["fasta", "fastq", "sam", "bam", "cram"];

    for input in &inputs {
        for output in &outputs {
            let in_format = input.split('.').next_back().unwrap();
            let mut cmd = Command::cargo_bin(BIN)?;
            let out_file = temp_dir
                .path()
                .join(format!("out_{}_{}.{}", in_format, output, output));

            cmd.args(vec![
                READS,
                input,
                "-n",
                "2",
                "-o",
                out_file.to_str().unwrap(),
                "-O",
                output,
            ]);

            // Verify that we also get a failure when writing to stdout (no -o flag)
            if (in_format == "fasta" || in_format == "fastq")
                && (output == &"sam" || output == &"bam" || output == &"cram")
            {
                let mut cmd_stdout = Command::cargo_bin(BIN)?;
                cmd_stdout.args(vec![READS, input, "-n", "2", "-O", output]);
                cmd_stdout.assert().failure();

                // FASTA/FASTQ cannot be converted to SAM/BAM/CRAM with -o flag
                cmd.assert().failure();
                continue;
            }

            cmd.assert().success();

            // Verify output format
            if output == &"fasta" {
                let content = fs::read_to_string(&out_file).unwrap();
                assert!(
                    content.starts_with(">"),
                    "Expected FASTA output for {} -> {}",
                    input,
                    output
                );
            } else if output == &"fastq" {
                let content = fs::read_to_string(&out_file).unwrap();
                assert!(
                    content.starts_with("@"),
                    "Expected FASTQ output for {} -> {}",
                    input,
                    output
                );
            } else if output == &"sam" {
                let content = fs::read_to_string(&out_file).unwrap();
                assert!(
                    content.starts_with("@HD")
                        || content.starts_with("@PG")
                        || content.contains("\t"),
                    "Expected SAM output for {} -> {}",
                    input,
                    output
                );
            } else if output == &"bam" || output == &"cram" {
                // Just check if it's not text
                let content = fs::read(&out_file).unwrap();
                assert!(
                    !content.starts_with(b">") && !content.starts_with(b"@\n"),
                    "Expected binary output for {} -> {}",
                    input,
                    output
                );
            }
        }
    }

    Ok(())
}

#[test]
fn one_pass_without_frac_raises_error() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/seed.fastq",
        "--one-pass",
        "-n",
        "5",
    ]);

    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("--one-pass requires --frac"));

    Ok(())
}

#[test]
fn one_pass_with_bases_raises_error() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/seed.fastq",
        "--one-pass",
        "-b",
        "100",
    ]);

    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("--one-pass requires --frac"));

    Ok(())
}

#[test]
fn one_pass_with_strict_raises_explanatory_error() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/seed.fastq",
        "--one-pass",
        "-f",
        "0.5",
        "--strict",
    ]);

    cmd.assert().failure().stderr(predicate::str::contains(
        "--one-pass cannot be combined with --strict",
    ));

    Ok(())
}

#[test]
fn one_pass_paired_read_count_mismatch_raises_error() -> Result<(), Box<dyn std::error::Error>> {
    // file1.fq.gz has 1 read, r2.fq.gz has 2 - a real count mismatch, detected mid-stream.
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/file1.fq.gz",
        "tests/cases/r2.fq.gz",
        "--one-pass",
        "-f",
        "0.5",
        "-o",
        "/tmp/one_pass_paired_mismatch_out1.fq",
        "-o",
        "/tmp/one_pass_paired_mismatch_out2.fq",
    ]);

    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("different numbers of reads"));

    Ok(())
}

#[test]
fn one_pass_with_bam_paired_input_raises_error() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/ubam/single_ubam.bam",
        "tests/cases/ubam/single_ubam.bam",
        "--one-pass",
        "-f",
        "0.5",
        "-o",
        "/tmp/one_pass_bam_paired_out1.fq",
        "-o",
        "/tmp/one_pass_bam_paired_out2.fq",
    ]);

    cmd.assert().failure().stderr(predicate::str::contains(
        "--one-pass only supports FASTA/FASTQ input",
    ));

    Ok(())
}

#[test]
fn one_pass_paired_keeps_mates_together_and_reports_once() -> Result<(), Box<dyn std::error::Error>>
{
    let temp_dir = tempfile::tempdir().unwrap();
    let out1 = temp_dir.path().join("out1.fastq");
    let out2 = temp_dir.path().join("out2.fastq");
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/r1.fq.gz",
        "tests/cases/r2.fq.gz",
        "--one-pass",
        "-f",
        "0.5",
        "-s",
        "1",
        "-o",
        out1.to_str().unwrap(),
        "-o",
        out2.to_str().unwrap(),
    ]);

    let stderr_output = cmd.output().unwrap().stderr;
    let stderr = String::from_utf8(stderr_output).unwrap();
    // reported once for the pair, not once per file
    assert_eq!(stderr.matches("reads seen").count(), 1);

    // r1/r2 use "/1" and "/2" mate-suffixed read names, so strip those before comparing - the
    // same *template* (read number) should be kept in both files, in the same order.
    let strip_mate_suffix = |content: String| -> Vec<String> {
        content
            .lines()
            .filter(|l| l.starts_with('@'))
            .map(|l| l.trim_end_matches("/1").trim_end_matches("/2").to_owned())
            .collect()
    };
    let ids1 = strip_mate_suffix(fs::read_to_string(&out1).unwrap());
    let ids2 = strip_mate_suffix(fs::read_to_string(&out2).unwrap());

    assert_eq!(ids1, ids2);
    assert!(!ids1.is_empty());

    Ok(())
}

#[test]
fn one_pass_with_bam_input_raises_error() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/ubam/single_ubam.bam",
        "--one-pass",
        "-f",
        "0.5",
    ]);

    cmd.assert().failure().stderr(predicate::str::contains(
        "--one-pass only supports FASTA/FASTQ input",
    ));

    Ok(())
}

#[test]
fn one_pass_keeping_zero_reads_warns_but_succeeds() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/seed.fastq",
        "--one-pass",
        "-f",
        "0",
        "-s",
        "1",
    ]);

    cmd.assert()
        .success()
        .stderr(predicate::str::contains("Kept 0 reads"));

    Ok(())
}

#[test]
fn one_pass_without_seed_logs_chosen_seed() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/seed.fastq",
        "--one-pass",
        "-f",
        "0.5",
    ]);

    cmd.assert()
        .success()
        .stderr(predicate::str::contains("Using seed:"));

    Ok(())
}

#[test]
fn one_pass_reports_kept_seen_and_fractions() -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempfile::tempdir().unwrap();
    let out = temp_dir.path().join("out.fastq");
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/seed.fastq",
        "--one-pass",
        "-f",
        "0.5",
        "-s",
        "1",
        "-o",
        out.to_str().unwrap(),
    ]);

    cmd.assert().success().stderr(
        predicate::str::contains("reads seen")
            .and(predicate::str::contains("requested fraction"))
            .and(predicate::str::contains("realised fraction")),
    );

    Ok(())
}

#[test]
fn one_pass_preserves_input_order_at_cli_level() -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempfile::tempdir().unwrap();
    let out = temp_dir.path().join("out.fastq");
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/seed.fastq",
        "--one-pass",
        "-f",
        "0.6",
        "-s",
        "7",
        "-o",
        out.to_str().unwrap(),
    ]);
    cmd.assert().success();

    let content = fs::read_to_string(&out).unwrap();
    let kept_ids: Vec<&str> = content.lines().filter(|l| l.starts_with('@')).collect();

    let input_content = fs::read_to_string("tests/cases/seed.fastq").unwrap();
    let input_order: Vec<&str> = input_content
        .lines()
        .filter(|l| l.starts_with('@'))
        .collect();
    // The kept ids, in the order they were written, should be exactly the subsequence of the
    // input's read order that was kept - not, say, a lexicographic sort of the ids.
    let expected_order: Vec<&str> = input_order
        .into_iter()
        .filter(|id| kept_ids.contains(id))
        .collect();

    assert_eq!(kept_ids, expected_order);
    assert!(!kept_ids.is_empty());

    Ok(())
}

#[test]
fn probability_shorthand_matches_frac_and_one_pass_output() -> Result<(), Box<dyn std::error::Error>>
{
    let temp_dir = tempfile::tempdir().unwrap();
    let long_form_out = temp_dir.path().join("long_form.fastq");
    let shorthand_out = temp_dir.path().join("shorthand.fastq");

    let mut long_form_cmd = Command::cargo_bin(BIN)?;
    long_form_cmd.args(vec![
        READS,
        "tests/cases/seed.fastq",
        "--one-pass",
        "-f",
        "0.5",
        "-s",
        "7",
        "-o",
        long_form_out.to_str().unwrap(),
    ]);
    long_form_cmd.assert().success();

    let mut shorthand_cmd = Command::cargo_bin(BIN)?;
    shorthand_cmd.args(vec![
        READS,
        "tests/cases/seed.fastq",
        "-p",
        "0.5",
        "-s",
        "7",
        "-o",
        shorthand_out.to_str().unwrap(),
    ]);
    shorthand_cmd.assert().success();

    let long_form_content = fs::read(&long_form_out).unwrap();
    let shorthand_content = fs::read(&shorthand_out).unwrap();
    assert_eq!(long_form_content, shorthand_content);

    Ok(())
}

#[test]
fn probability_shorthand_matches_frac_and_one_pass_output_for_paired_input(
) -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempfile::tempdir().unwrap();
    let long_form_out1 = temp_dir.path().join("long_form_r1.fastq");
    let long_form_out2 = temp_dir.path().join("long_form_r2.fastq");
    let shorthand_out1 = temp_dir.path().join("shorthand_r1.fastq");
    let shorthand_out2 = temp_dir.path().join("shorthand_r2.fastq");

    let mut long_form_cmd = Command::cargo_bin(BIN)?;
    long_form_cmd.args(vec![
        READS,
        "tests/cases/seed_r1.fastq",
        "tests/cases/seed_r2.fastq",
        "--one-pass",
        "-f",
        "0.5",
        "-s",
        "7",
        "-o",
        long_form_out1.to_str().unwrap(),
        "-o",
        long_form_out2.to_str().unwrap(),
    ]);
    long_form_cmd.assert().success();

    let mut shorthand_cmd = Command::cargo_bin(BIN)?;
    shorthand_cmd.args(vec![
        READS,
        "tests/cases/seed_r1.fastq",
        "tests/cases/seed_r2.fastq",
        "-p",
        "0.5",
        "-s",
        "7",
        "-o",
        shorthand_out1.to_str().unwrap(),
        "-o",
        shorthand_out2.to_str().unwrap(),
    ]);
    shorthand_cmd.assert().success();

    assert_eq!(
        fs::read(&long_form_out1).unwrap(),
        fs::read(&shorthand_out1).unwrap()
    );
    assert_eq!(
        fs::read(&long_form_out2).unwrap(),
        fs::read(&shorthand_out2).unwrap()
    );

    Ok(())
}

#[test]
fn probability_shorthand_interprets_values_over_one_as_percentages(
) -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempfile::tempdir().unwrap();
    let percentage_out = temp_dir.path().join("percentage.fastq");
    let fraction_out = temp_dir.path().join("fraction.fastq");

    let mut percentage_cmd = Command::cargo_bin(BIN)?;
    percentage_cmd.args(vec![
        READS,
        "tests/cases/seed.fastq",
        "-p",
        "50",
        "-s",
        "7",
        "-o",
        percentage_out.to_str().unwrap(),
    ]);
    percentage_cmd.assert().success();

    let mut fraction_cmd = Command::cargo_bin(BIN)?;
    fraction_cmd.args(vec![
        READS,
        "tests/cases/seed.fastq",
        "-p",
        "0.5",
        "-s",
        "7",
        "-o",
        fraction_out.to_str().unwrap(),
    ]);
    fraction_cmd.assert().success();

    let percentage_content = fs::read(&percentage_out).unwrap();
    let fraction_content = fs::read(&fraction_out).unwrap();
    assert_eq!(percentage_content, fraction_content);

    Ok(())
}

#[test]
fn probability_and_frac_not_allowed() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/seed.fastq",
        "-p",
        "0.5",
        "-f",
        "0.5",
    ]);

    cmd.assert().failure();

    Ok(())
}

#[test]
fn probability_and_num_not_allowed() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/seed.fastq",
        "-p",
        "0.5",
        "-n",
        "5",
    ]);

    cmd.assert().failure();

    Ok(())
}

#[test]
fn probability_and_bases_not_allowed() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/seed.fastq",
        "-p",
        "0.5",
        "-b",
        "100",
    ]);

    cmd.assert().failure();

    Ok(())
}

#[test]
fn probability_and_genome_size_not_allowed() -> Result<(), Box<dyn std::error::Error>> {
    // `-c` is also passed so this fails specifically because `-g` conflicts with `-p`, not
    // because `-g` is missing its required `-c` companion.
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/seed.fastq",
        "-p",
        "0.5",
        "-g",
        "5m",
        "-c",
        "10",
    ]);

    cmd.assert().failure();

    Ok(())
}

#[test]
fn probability_and_coverage_not_allowed() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/seed.fastq",
        "-p",
        "0.5",
        "-c",
        "10",
        "-g",
        "5m",
    ]);

    cmd.assert().failure();

    Ok(())
}

#[test]
fn probability_and_strict_not_allowed() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        READS,
        "tests/cases/seed.fastq",
        "-p",
        "0.5",
        "--strict",
    ]);

    cmd.assert().failure();

    Ok(())
}

#[test]
fn probability_alone_does_not_require_genome_size_or_coverage(
) -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![READS, "tests/cases/seed.fastq", "-p", "0.5"]);

    cmd.assert().success();

    Ok(())
}

#[test]
fn probability_help_text_documents_the_equivalence_and_approximation(
) -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![READS, "--help"]);

    cmd.assert().success().stdout(
        predicate::str::contains("--frac <FLOAT> --one-pass")
            .and(predicate::str::contains("is approximate")),
    );

    Ok(())
}
