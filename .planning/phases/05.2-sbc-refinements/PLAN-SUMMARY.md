# Plan: SBC Refinements and Documentation

Final refinements for SBC support: fixed BAM-to-FASTA bug, updated CLI help, and documentation.

## Accomplishments
- Fixed bug where BAM to FASTA conversion would output FASTQ format.
- Added `is_fasta_path` helper in `src/reads.rs` for better format inference.
- Updated CLI help strings in `src/reads.rs` to emphasize unaligned SAM/BAM/CRAM.
- Updated `README.md` and planning docs to specify "unaligned" SBC support.
- Added `reads_bam_to_fasta` integration test in `tests/main.rs`.

## Verification Results
- `cargo test reads_bam_to_fasta` PASSED.
- All integration and unit tests passing.
