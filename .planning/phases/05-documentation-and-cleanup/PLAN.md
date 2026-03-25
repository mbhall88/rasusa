---
phase: 05-documentation-and-cleanup
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - README.md
  - src/source.rs
  - src/reads.rs
  - src/fastx.rs
  - .planning/STATE.md
  - .planning/ROADMAP.md
autonomous: true
requirements:
  - REQ-DOC-README
  - REQ-CLEANUP
  - REQ-PROJECT-COMPLETE
must_haves:
  truths:
    - "README.md explains SAM/BAM/CRAM support for the reads subcommand"
    - "README.md specifies the new extension-based format inference for outputs in reads"
    - "Performance check on a BAM file runs successfully"
    - "STATE.md and ROADMAP.md reflect project completion"
  artifacts:
    - path: "README.md"
      provides: "Updated documentation for SAM/BAM/CRAM support in reads"
    - path: ".planning/STATE.md"
      provides: "Final project state"
    - path: ".planning/ROADMAP.md"
      provides: "Completed roadmap"
  key_links:
    - from: "README.md"
      to: "reads subcommand"
      via: "Description of new SAM/BAM/CRAM support"
---

<objective>
Update project documentation to reflect the new SAM/BAM/CRAM support in the `reads` subcommand, perform a final code review and performance sanity check, and finalize the project state to indicate full completion.

Purpose: Ensure the new functionality is well-documented and the codebase is production-ready.
Output: Updated README.md, reviewed code, and finalized project state files.
</objective>

<execution_context>
@$HOME/.gemini/get-shit-done/workflows/execute-plan.md
@$HOME/.gemini/get-shit-done/templates/summary.md
</execution_context>

<context>
@.planning/ROADMAP.md
@.planning/STATE.md
@README.md
@src/source.rs
@src/reads.rs
@src/fastx.rs
</context>

<tasks>

<task type="auto">
  <name>Task 1: Update README.md for SAM/BAM/CRAM support in reads</name>
  <files>README.md</files>
  <action>
    - Update "Basic usage - reads" to mention SAM/BAM/CRAM support and add an example: `rasusa reads --coverage 30 --genome-size 4.6mb in.bam`.
    - Update "Required parameters -> Input" to reflect that `reads` now accepts SAM, BAM, and CRAM files in addition to FASTA/FASTQ.
    - Update "Optional parameters -> Output" for `reads` to explain that output format is inferred from the extension (e.g., .bam, .sam, .cram) and that conversion is supported.
    - Update "Optional parameters -> Output compression/format" for `reads` to include alignment format options (b/bam, c/cram, s/sam) and mention that these override extension inference.
  </action>
  <verify>
    <automated>grep -A 10 "Basic usage - reads" README.md | grep ".bam"</automated>
  </verify>
  <done>README.md accurately documents the new SAM/BAM/CRAM support for the reads subcommand across all relevant sections.</done>
</task>

<task type="auto">
  <name>Task 2: Final code review and performance check</name>
  <files>src/source.rs, src/reads.rs, src/fastx.rs</files>
  <action>
    - Perform a final code review of `src/source.rs`, `src/reads.rs`, and `src/fastx.rs` to ensure consistency, proper error handling, and clean implementation of the `RecordSource` abstraction.
    - Run a performance sanity check by subsampling a BAM file (e.g., `tests/cases/test.bam`) with the `reads` subcommand and ensuring it completes efficiently.
    - Verify that all tests pass: `cargo test`.
  </action>
  <verify>
    <automated>cargo test && cargo run -- reads --coverage 10 --genome-size 100k tests/cases/test.bam > /dev/null</automated>
  </verify>
  <done>Code is reviewed and polished, and the performance sanity check on BAM input confirms the implementation is sound.</done>
</task>

<task type="auto">
  <name>Task 3: Finalize project state and roadmap</name>
  <files>.planning/STATE.md, .planning/ROADMAP.md</files>
  <action>
    - Update `.planning/STATE.md` to mark Phase 5 as completed and indicate the project is now in its final state.
    - Update `.planning/ROADMAP.md` to mark Phase 5 as [x] and ensure all previous phases are correctly marked as completed.
  </action>
  <verify>
    <automated>grep "\[x\] Phase 5" .planning/ROADMAP.md</automated>
  </verify>
  <done>Project state and roadmap are fully up to date, reflecting 100% completion of the SAM/BAM/CRAM support project.</done>
</task>

</tasks>

<verification>
- README.md has been updated with the new documentation.
- Final code review completed for core modules.
- Performance sanity check on BAM file passes.
- STATE.md and ROADMAP.md show the project as completed.
</verification>

<success_criteria>
- Clear and accurate documentation for users to utilize the new SAM/BAM/CRAM support in the `reads` subcommand.
- Verified and clean implementation of the core logic.
- Project management files correctly reflect the finished state of the project.
</success_criteria>

<output>
After completion, create `.planning/phases/05-documentation-and-cleanup/05-01-SUMMARY.md`
</output>
