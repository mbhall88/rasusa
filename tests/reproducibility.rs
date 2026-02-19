use assert_cmd::Command;

const BIN: &str = "rasusa";

#[test]
fn reproducibility_reads_seeds_1_to_5() -> Result<(), Box<dyn std::error::Error>> {
    let cases = vec![
        (
            1,
            vec![
                "@read1", "@read2", "@read6", "@read7", "@read9", "@read10", "@read11", "@read12",
                "@read13", "@read14",
            ],
        ),
        (
            2,
            vec![
                "@read3", "@read4", "@read5", "@read6", "@read7", "@read9", "@read10", "@read11",
                "@read12", "@read15",
            ],
        ),
        (
            3,
            vec![
                "@read1", "@read2", "@read3", "@read4", "@read8", "@read10", "@read11", "@read13",
                "@read14", "@read15",
            ],
        ),
        (
            4,
            vec![
                "@read1", "@read2", "@read3", "@read4", "@read5", "@read6", "@read11", "@read13",
                "@read15", "@read16",
            ],
        ),
        (
            5,
            vec![
                "@read2", "@read5", "@read7", "@read8", "@read9", "@read10", "@read12", "@read13",
                "@read14", "@read15",
            ],
        ),
    ];

    for (seed, expected_reads) in cases {
        let mut cmd = Command::cargo_bin(BIN)?;
        cmd.args(vec![
            "reads",
            "tests/cases/seed.fastq",
            "-s",
            &seed.to_string(),
            "-n",
            "10",
        ]);

        let output = cmd.output()?;
        let stdout = String::from_utf8(output.stdout)?;
        let actual_reads: Vec<&str> = stdout.lines().filter(|l| l.starts_with('@')).collect();

        assert_eq!(
            actual_reads, expected_reads,
            "Reproduction failed for seed {}",
            seed
        );
    }

    Ok(())
}

#[test]
fn reproducibility_aln_seed_42() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![
        "aln",
        "tests/cases/test.bam",
        "-s",
        "42",
        "-c",
        "5",
        "-o",
        "test_aln_out.sam",
    ]);

    cmd.assert().success();

    let output = std::fs::read_to_string("test_aln_out.sam")?;
    let mut actual_reads: Vec<String> = output
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(|l| l.split('\t').next().unwrap().to_string())
        .collect();
    actual_reads.sort();
    actual_reads.dedup();

    let expected_output = std::fs::read_to_string("tests/cases/baseline_aln_42.txt")?;
    let expected_reads: Vec<String> = expected_output.lines().map(|l| l.to_string()).collect();

    assert_eq!(actual_reads, expected_reads);

    // Clean up
    let _ = std::fs::remove_file("test_aln_out.sam");

    Ok(())
}
