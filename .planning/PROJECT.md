# Project: rasusa

**Description:** Randomly subsample reads or alignments.
**Repository:** https://github.com/mbhall88/rasusa
**Primary Language:** Rust (v1.84.1)

## High-Level Goals
- Efficiently subsample biological sequence files (FASTA/FASTQ/unaligned SAM/BAM/CRAM).
- Subsample aligned SAM/BAM/CRAM files while maintaining per-position coverage depth.
- Provide a robust, user-friendly CLI.

## Current State
- Version 3.0.0.
- Established support for multiple formats (FASTX and unaligned SAM/BAM/CRAM) via `noodles` and `needletail`.
- Known technical debt in `src/alignment.rs` and double-pass I/O in `src/reads.rs`.

## Strategy
- Implement new features with a focus on efficiency and performance.
- Refactor fragile areas as they are touched by new functionality.
- Maintain high test coverage using `assert_cmd` and `predicates`.
- **Note:** Use the `SEGMENTED` flag (bit 0x1) in SAM/BAM/CRAM records to detect paired-end data.
