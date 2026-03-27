---
phase: 05.1-fix-ubam-output
plan: 01
type: execute
wave: 1
depends_on: []
files_modified: [src/reads.rs, src/alignment.rs, src/main.rs]
autonomous: true
requirements: [UBAM-MATCH-FORMAT, UBAM-CONVERT-FASTQ, UBAM-HEADER-PG, UBAM-VALIDATE-UNMAPPED]

must_haves:
  truths:
    - "Output format matches input SAM/BAM/CRAM by default"
    - "Conversion from SBC to FASTQ works via extension"
    - "SAM/BAM/CRAM headers contain 'rasusa' PG entry"
    - "Mapped reads in SBC input trigger specific error message"
  artifacts:
    - path: "src/alignment.rs"
      provides: "SBC validation and header manipulation"
  key_links:
    - from: "rasusa reads"
      to: "SBC output"
      via: "output extension detection"
---

<objective>
Refine unaligned SAM/BAM/CRAM (SBC) support in the `reads` subcommand to ensure format parity, header traceability, and data integrity.
</objective>

<execution_context>
@$HOME/.gemini/get-shit-done/workflows/execute-plan.md
</execution_context>

<context>
@.planning/PROJECT.md
@.planning/ROADMAP.md
@src/alignment.rs
@src/reads.rs
</context>

<tasks>

<task type="auto">
  <name>Task 1: Refine SBC Output and Headers</name>
  <files>src/alignment.rs, src/reads.rs</files>
  <action>
    - Update the writer logic to default to the input SBC format if no output format is explicitly specified.
    - Implement logic to allow SBC to FASTQ conversion based on output file extension.
    - Explicitly block FASTQ to SBC conversion if attempted.
    - Add a 'rasusa' @PG (Program) entry to the SAM/BAM/CRAM header, including version information, following standard alignment tool conventions.
  </action>
  <verify>
    <automated>cargo test --test test_writer</automated>
  </verify>
  <done>SBC outputs have correct headers and formats match expectations.</done>
</task>

<task type="auto">
  <name>Task 2: Implement Unmapped Validation</name>
  <files>src/alignment.rs, src/main.rs</files>
  <action>
    - Add a check for each record in SBC input to verify the unmapped flag (0x04) is set.
    - If a mapped read is detected, raise a fatal error with the message: "Error: Mapped read detected, please use `rasusa aln` for aligned data".
  </action>
  <verify>
    <automated>cargo run -- reads --input tests/cases/test.bam --coverage 10 --output out.bam || true</automated>
  </verify>
  <done>Mapped reads correctly trigger the specified error message.</done>
</task>

<task type="auto">
  <name>Task 3: Integration Testing with uBAM cases</name>
  <files>tests/main.rs</files>
  <action>
    - Add integration tests using the test cases in `tests/cases/ubam/` (paired_interleave_ubam.bam, single_ubam.bam, etc.).
    - Verify that `rasusa reads` processes these files correctly and preserves (or converts) formats as requested.
  </action>
  <verify>
    <automated>cargo test --test main ubam</automated>
  </verify>
  <done>All uBAM test cases pass with expected output formats and headers.</done>
</task>

</tasks>

<success_criteria>
- All records in `reads` SBC output are verified unmapped.
- SBC headers contain rasusa PG entry.
- Format parity is maintained by default.
- Integration tests for uBAM cases pass.
</success_criteria>
