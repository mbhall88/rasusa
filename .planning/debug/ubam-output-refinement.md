---
status: investigating
trigger: "Investigate and refine unaligned SAM/BAM/CRAM (ubam) support in the `reads` subcommand."
created: 2024-05-24T12:00:00Z
updated: 2024-05-24T12:00:00Z
---

## Current Focus

hypothesis: The current `reads` subcommand implementation defaults to FASTX output because it uses a FASTX-specific writer or doesn't detect SAM/BAM/CRAM input formats for output format matching.
test: Examine `src/cli.rs` and `src/main.rs` to see how the `reads` subcommand is handled.
expecting: Find logic that forces FASTX output regardless of input format.
next_action: gather initial evidence by reading relevant source files.

## Symptoms

expected: 
- SAM/BAM/CRAM input -> matching SAM/BAM/CRAM output by default.
- Conversion from SAM/BAM/CRAM to FASTQ allowed if specified by output extension.
- SAM/BAM/CRAM output headers include a 'rasusa' program entry (Standard aln entry style).
- Error "Error: Mapped read detected, please use `rasusa aln` for aligned data" if flag 0x04 (UNMAPPED) is not set.
actual: 
- Output is always FASTA/FASTQ.
- No program entry added.
- No verification of the unmapped flag.
errors: N/A (Missing feature/logic error)
reproduction: Use test cases in `tests/cases/ubam/`.
started: Refinement of Phase 1 functionality.

## Eliminated

## Evidence

## Resolution

root_cause: 
fix: 
verification: 
files_changed: []
