---
phase: 05.2-sbc-refinements
plan: 02
type: execute
wave: 1
depends_on: []
files_modified: [src/cli.rs, src/reads.rs, src/alignment.rs, tests/main.rs]
autonomous: true
requirements: [CLI-FLAG-REFACTOR]

must_haves:
  truths:
    - "`--compress-type` controls output compression instead of `--output-type`"
    - "New `--output-format` explicit flag allows overriding format inference"
    - "Flags correctly use `-Z` for compression and `-O` for format"
  artifacts:
    - path: "src/cli.rs"
      provides: "OutputFormat enum definition"
    - path: "src/reads.rs"
      provides: "Reads struct with updated flags"
    - path: "src/alignment.rs"
      provides: "Alignment struct with updated flags"
---

<objective>
Refactor CLI flags to explicitly separate file format from compression type.
</objective>

<execution_context>
@$HOME/.gemini/get-shit-done/workflows/execute-plan.md
</execution_context>

<context>
@.planning/phases/05.2-sbc-refinements/05.2-CONTEXT.md
@src/cli.rs
@src/reads.rs
@src/alignment.rs
</context>

<tasks>

<task type="auto">
  <name>Task 1: Add OutputFormat Enum</name>
  <files>src/cli.rs</files>
  <action>
    - Define a new enum `OutputFormat` with variants: `Fasta`, `Fastq`, `Sam`, `Bam`, `Cram`.
    - Derive `clap::ValueEnum`, `Clone`, `Debug` for `OutputFormat`.
  </action>
  <read_first>src/cli.rs</read_first>
  <acceptance_criteria>
    - `src/cli.rs` contains `pub enum OutputFormat` with the specified variants.
  </acceptance_criteria>
  <done>OutputFormat enum is available for clap parsing.</done>
</task>

<task type="auto">
  <name>Task 2: Refactor Reads Struct Flags</name>
  <files>src/reads.rs</files>
  <action>
    - Rename `pub output_type: Option<niffler::compression::Format>` to `pub compress_type: Option<niffler::compression::Format>`.
    - Change its `clap` attribute to use `short = 'Z'`, `long = "compress-type"`.
    - Add `#[clap(short = 'O', long = "output-format", value_enum)] pub output_format: Option<crate::cli::OutputFormat>,`.
    - Update all usages of `self.output_type` in `src/reads.rs` to `self.compress_type`.
  </action>
  <read_first>src/reads.rs</read_first>
  <acceptance_criteria>
    - `grep -q "pub compress_type: Option<niffler::compression::Format>" src/reads.rs` succeeds.
    - `grep -q "pub output_format: Option<crate::cli::OutputFormat>" src/reads.rs` succeeds.
    - `cargo check` compiles successfully (or fails only on `source.rs` which is fixed in PLAN 01).
  </acceptance_criteria>
  <done>Reads subcommand flags are refactored.</done>
</task>

<task type="auto">
  <name>Task 3: Refactor Alignment Struct Flags</name>
  <files>src/alignment.rs</files>
  <action>
    - Rename `pub output_type: Option<Format>` to `pub output_format: Option<Format>`.
    - Update the `clap` attribute to keep `short='O', long="output-format"`.
    - Update all usages of `self.output_type` to `self.output_format` in `src/alignment.rs`.
  </action>
  <read_first>src/alignment.rs</read_first>
  <acceptance_criteria>
    - `grep -q "pub output_format: Option<Format>" src/alignment.rs` succeeds.
    - `cargo check` compiles successfully for alignment.
  </acceptance_criteria>
  <done>Alignment subcommand flags are refactored consistently.</done>
</task>

</tasks>

<success_criteria>
- `OutputFormat` enum exists in `src/cli.rs`.
- `Reads` subcommand uses `--compress-type` (`-Z`) and `--output-format` (`-O`).
- `Alignment` subcommand uses `--output-format` (`-O`).
- `cargo build` succeeds after these changes.
</success_criteria>
