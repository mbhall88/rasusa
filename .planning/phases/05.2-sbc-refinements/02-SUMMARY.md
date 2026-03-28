# Plan 02: CLI Flag Refinement for SBC

Refactored CLI flags to explicitly separate file format from compression type.

## Accomplishments
- Added `OutputFormat` enum to `src/cli.rs`.
- Renamed `--output-type` to `--compress-type` (short flag `-Z`) in `src/reads.rs`.
- Added `--output-format` (short flag `-O`) to `src/reads.rs`.
- Renamed `output_type` to `output_format` in `src/alignment.rs` for consistency.
- Updated `filter_reads_into` signature to take an `is_fasta` boolean to control format output explicitly.

## Verification Results
- Build passes (`cargo check`).
- Tests in `src/source.rs` and `src/fastx.rs` updated and passing.
