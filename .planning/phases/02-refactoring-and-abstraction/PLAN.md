---
phase: 02-refactoring-and-abstraction
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - src/source.rs
  - src/fastx.rs
  - src/reads.rs
  - src/main.rs
autonomous: true
requirements:
  - REQ-REFACTOR
must_haves:
  truths:
    - "System routes input files to the correct format handler dynamically"
    - "Fastx object only handles FASTA/FASTQ logic"
    - "AlignmentSource object handles SAM/BAM/CRAM logic"
    - "Output writer is decoupled from input source"
  artifacts:
    - path: "src/source.rs"
      provides: "RecordSource trait and determine_record_source factory"
    - path: "src/fastx.rs"
      provides: "Fastx struct implementing RecordSource"
  key_links:
    - from: "src/reads.rs"
      to: "src/source.rs"
      via: "Box<dyn RecordSource> dynamic dispatch"
---

<objective>
Refactor `src/reads.rs` and `src/fastx.rs` to use a generic `RecordSource` abstraction.

Purpose: To resolve the Single Responsibility Principle violation by removing format-specific branching logic (SAM/BAM/CRAM vs FASTA/FASTQ) from `Fastx` and standardizing execution in `src/reads.rs` using dynamic dispatch (`Box<dyn RecordSource>`).
Output: A new `src/source.rs` module containing the `RecordSource` trait and factory function, `Fastx` and `AlignmentSource` implementations, and a refactored `reads` subcommand execution flow.
</objective>

<context>
@.planning/ROADMAP.md
@.planning/phases/02-refactoring-and-abstraction/RESEARCH.md
@src/reads.rs
@src/fastx.rs
</context>

<tasks>

<task type="auto">
  <name>Task 1: Define RecordSource trait and standalone writer</name>
  <files>src/source.rs, src/fastx.rs, src/main.rs</files>
  <action>
    - Create `src/source.rs` and declare `pub trait RecordSource`. Make methods object-safe: `fn read_lengths(&self) -> Result<Vec<u32>, crate::fastx::FastxError>` and `fn filter_reads_into(&self, reads_to_keep: &[bool], nb_reads_keep: usize, write_to: &mut dyn std::io::Write, output_format: Option<noodles_util::alignment::io::Format>) -> Result<usize, crate::fastx::FastxError>`.
    - Extract output creation logic (`Fastx::create`) from `src/fastx.rs` into a standalone function `pub fn create_output_writer(...) -> Result<Box<dyn std::io::Write>, FastxError>` in `src/fastx.rs` (or `src/source.rs`). Remove the `create` method from `Fastx`.
    - Update unit tests in `src/fastx.rs` (like `create_valid_output_file_and_can_write_to_it` and others testing file creation) to use the new standalone `create_output_writer` function instead of `Fastx::create`.
    - Register `mod source;` in `src/main.rs`.
  </action>
  <verify>
    <automated>cargo check</automated>
  </verify>
  <done>RecordSource trait is defined with object-safe methods, output writer creation is extracted, and creation-related tests are updated.</done>
</task>

<task type="auto">
  <name>Task 2: Implement RecordSource for format handlers</name>
  <files>src/fastx.rs, src/source.rs</files>
  <action>
    - In `src/fastx.rs`: Implement `RecordSource` for `Fastx`. Remove any SAM/BAM/CRAM branching logic (`match ext.as_str()`) from `Fastx` methods. It should only handle FASTA/FASTQ reading and writing.
    - In `src/source.rs` (or `src/alignment.rs` if appropriate): Create an `AlignmentSource` struct (holding the file path). Implement `RecordSource` for `AlignmentSource`, moving the SAM/BAM/CRAM parsing and filtering logic from the old `Fastx` implementation into this struct.
    - Implement a factory function `pub fn determine_record_source(path: &std::path::Path) -> Box<dyn RecordSource>` in `src/source.rs` that returns either `Box<AlignmentSource>` or `Box<Fastx>` based on the file extension.
    - Migrate SAM/BAM/CRAM specific unit tests from `src/fastx.rs` (e.g., `get_read_lengths_for_bam`) to the new module containing `AlignmentSource` so they test the new implementation instead of `Fastx`. Update them to use `AlignmentSource` instead of `Fastx`.
  </action>
  <verify>
    <automated>cargo check</automated>
  </verify>
  <done>Fastx and AlignmentSource cleanly implement RecordSource without cross-format logic branching, and BAM/SAM tests are correctly migrated.</done>
</task>

<task type="auto">
  <name>Task 3: Refactor reads.rs execution flow</name>
  <files>src/reads.rs</files>
  <action>
    - Update `src/reads.rs` to dynamically detect format and use `Box<dyn RecordSource>` instead of hardcoding the `Fastx` struct.
    - Call the new `determine_record_source(&self.input[0])` factory for input files instead of `Fastx::from_path`.
    - Replace the call to `Fastx::create` with the new standalone `create_output_writer` function to set up `output_handle`.
    - Ensure `read_lengths()` and `filter_reads_into()` are invoked on the trait object (`input_source.read_lengths()`, etc.).
  </action>
  <verify>
    <automated>cargo test</automated>
  </verify>
  <done>Reads subcommand builds successfully, uses dynamic dispatch, and all existing unit/integration tests pass.</done>
</task>

</tasks>

<verification>
- `cargo check` passes without warnings about trait object safety.
- `cargo test` passes, indicating no regressions in FASTQ/FASTA support.
</verification>

<success_criteria>
- `Fastx` struct no longer contains noodles/alignment logic.
- `src/reads.rs` utilizes `Box<dyn RecordSource>` instead of `Fastx` type directly for inputs.
- Standalone `create_output_writer` cleanly handles output file stream creation.
</success_criteria>

<output>
After completion, create `.planning/phases/02-refactoring-and-abstraction/02-01-SUMMARY.md`
</output>
