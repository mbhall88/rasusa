---
phase: 04-validation-and-testing
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - src/source.rs
  - tests/main.rs
autonomous: true
requirements:
  - REQ-TEST-UNIT
  - REQ-TEST-PAIRED
  - REQ-FORMAT-CONV
  - REQ-TEST-INT
  - REQ-SEED
must_haves:
  truths:
    - "AlignmentSource::read_lengths correctly sums lengths for paired-end segments"
    - "AlignmentSource::filter_reads_into maintains paired-end integrity"
    - "AlignmentSource can convert between SAM, BAM, CRAM, and FASTQ/FASTA formats"
    - "rasusa reads subcommand correctly identifies and processes SBC input files"
    - "rasusa reads with fixed seed produces bit-identical results across multiple runs"
  artifacts:
    - path: "src/source.rs"
      provides: "Comprehensive unit tests for AlignmentSource"
    - path: "tests/main.rs"
      provides: "New integration tests for SBC and reproducibility"
  key_links:
    - from: "src/source.rs"
      to: "noodles"
      via: "Alignment records parsing and writing"
    - from: "rasusa cli"
      to: "subsampler"
      via: "Fixed seed for reproducibility"
---

<objective>
Comprehensive testing and validation of SAM/BAM/CRAM support in the `reads` subcommand, including unit tests, integration tests, format conversion, and reproducibility verification.

Purpose: To ensure the robustness and accuracy of the new format support and maintain high reliability for genomic data subsampling.
Output: Updated `src/source.rs` with unit tests and updated `tests/main.rs` with integration tests.
</objective>

<execution_context>
@$HOME/.gemini/get-shit-done/workflows/execute-plan.md
@$HOME/.gemini/get-shit-done/templates/summary.md
</execution_context>

<context>
@.planning/ROADMAP.md
@src/source.rs
@tests/main.rs
@tests/cases/test.bam
@tests/cases/test.paired.bam
</context>

<tasks>

<task type="auto">
  <name>Task 1: Unit tests for AlignmentSource (reading, filtering, conversion)</name>
  <files>src/source.rs</files>
  <action>
    - Add unit tests to `src/source.rs` to verify `read_lengths` and `filter_reads_into` for all alignment formats (SAM, BAM, CRAM).
    - Specifically verify paired-end integrity for single-file paired BAM (summing lengths and keeping/discarding pairs consistently).
    - Implement unit tests for format conversions: SBC to SBC (e.g., BAM to SAM) and SBC to FASTX (BAM to FASTQ/FASTA).
    - Verify that output records include correct names, sequences, and quality scores when applicable.
  </action>
  <verify>
    <automated>cargo test source::tests</automated>
  </verify>
  <done>Core alignment logic is fully verified at the unit level, covering reading, filtering, and format conversion.</done>
</task>

<task type="auto">
  <name>Task 2: Integration tests for reads subcommand and reproducibility</name>
  <files>tests/main.rs</files>
  <action>
    - Extend `tests/main.rs` with integration tests for `rasusa reads` using SAM, BAM, and CRAM as input.
    - Test auto-detection and manual override (`-O`) of output formats.
    - Implement reproducibility tests that run `rasusa reads` with a fixed `--seed` multiple times and compare outputs for identical content (checksums or bitwise comparison).
    - Verify that reproducibility holds for both FASTX and SBC input formats.
  </action>
  <verify>
    <automated>cargo test tests::main</automated>
  </verify>
  <done>Integration tests confirm end-to-end functionality and bit-identical reproducibility for all supported formats.</done>
</task>

</tasks>

<verification>
- `cargo test` passes all tests (unit and integration).
</verification>

<success_criteria>
- All segments of a pair are consistently kept or discarded in subsampling.
- BAM/SAM/CRAM files can be successfully processed and converted.
- Fixed seed ensures identical output on repeated runs.
- `rasusa reads` CLI behavior matches expectations for all formats.
</success_criteria>

<output>
After completion, create `.planning/phases/04-validation-and-testing/04-01-SUMMARY.md`
</output>
