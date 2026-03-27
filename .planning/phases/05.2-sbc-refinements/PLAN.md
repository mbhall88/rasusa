---
phase: 05.2-sbc-refinements
plan: 01
type: execute
wave: 1
depends_on: []
files_modified: [src/reads.rs, src/source.rs, README.md, .planning/PROJECT.md, .planning/REQUIREMENTS.md]
autonomous: true
requirements: [FIX-BAM-TO-FASTA, CLI-HELP-REFINEMENT, DOC-UNALIGNED-EMPHASIS]

must_haves:
  truths:
    - "BAM to .fasta output is valid FASTA (no quality scores)"
    - "CLI help mentions unaligned SAM/BAM/CRAM for reads subcommand"
    - "Documentation emphasizes 'unaligned' for SBC support"
  artifacts:
    - path: "src/source.rs"
      provides: "RecordSource implementation for SBC"
    - path: "src/reads.rs"
      provides: "CLI definition and runner logic"
---

<objective>
Final refinements for SBC support: fix BAM-to-FASTA bug, update CLI help, and emphasize "unaligned" SBC support in documentation.
</objective>

<execution_context>
@$HOME/.gemini/get-shit-done/workflows/execute-plan.md
</execution_context>

<context>
@.planning/PROJECT.md
@.planning/ROADMAP.md
@src/source.rs
@src/reads.rs
</context>

<tasks>

<task type="auto">
  <name>Task 1: Fix BAM to FASTA conversion</name>
  <files>src/source.rs</files>
  <action>
    - Investigate why `AlignmentSource::filter_reads_into` might be outputting FASTQ when FASTA is expected.
    - Specifically, check the extension-based inference logic in `src/reads.rs` and how it's passed to `filter_reads_into`.
    - Fix: Ensure that if `output_format` is `None` (meaning FASTX), we check the output path extension to decide between FASTA and FASTQ, OR ensure `noodles` records are written as FASTA if the extension is `.fasta`/`.fa`.
    - NOTE: The user says "When the input is BAM and I set -o test.fasta the output is actually fastq".
  </action>
  <verify>
    <automated>Add a new integration test in `tests/main.rs` called `reads_bam_to_fasta` that verifies the output starts with `>` and does not contain `+` or quality lines.</automated>
  </verify>
  <done>BAM to FASTA conversion produces valid FASTA content.</done>
</task>

<task type="auto">
  <name>Task 2: Refine CLI Options and Help</name>
  <files>src/reads.rs</files>
  <action>
    - Update `Reads` struct field comments to explicitly mention "unaligned SAM/BAM/CRAM".
    - Update `--output-type` help to mention it's for FASTA/FASTQ compression.
    - Clarify that `-o` extension is used for alignment format inference (SAM/BAM/CRAM).
    - Ensure help menu is clear and professional.
  </action>
  <verify>
    <automated>cargo run -- reads --help</automated>
  </verify>
  <done>CLI help is updated and accurate.</done>
</task>

<task type="auto">
  <name>Task 3: Documentation Update (Unaligned Emphasis)</name>
  <files>README.md, .planning/PROJECT.md, .planning/REQUIREMENTS.md, .planning/ROADMAP.md</files>
  <action>
    - Search for all occurrences of "SAM/BAM/CRAM" or "alignment files" in the context of the `reads` subcommand.
    - Prepend/qualify with "unaligned" where missing.
    - Update README usage examples to use "unaligned" terminology.
  </action>
  <verify>
    <automated>grep -i "unaligned" README.md</automated>
  </verify>
  <done>All documentation consistently emphasizes unaligned SBC support.</done>
</task>

</tasks>

<success_criteria>
- BAM to FASTA conversion works correctly.
- CLI help explicitly documents unaligned SBC output support.
- Documentation consistently refers to "unaligned" SAM/BAM/CRAM.
- All integration tests pass.
</success_criteria>
