# Plan 01-01 Summary: Core SBC reading and length gathering

## Tasks Completed
1. **Implement SBC read lengths gathering:** Modified `src/fastx.rs` to detect SAM/BAM/CRAM file extensions and utilize `noodles_util::alignment::io::reader::Builder` to read the sequence length of each record. Added `FastxError::AlignmentReadError` for error handling and added a new unit test `get_read_lengths_for_bam`.
2. **Ensure correct target base calculation:** Verified that by extending `Fastx::read_lengths()`, the `target_total_bases` calculation in `src/reads.rs` naturally inherits the correct values for SBC files since it relies on the returned vector of read lengths.

## Verification
- Code compiles successfully.
- `cargo test --bin rasusa fastx::tests` passes all 15 tests, including the new `get_read_lengths_for_bam` test which confirms it can successfully read lengths from a BAM file.

## Conclusion
The `reads` subcommand can now parse SAM, BAM, and CRAM files and accurately gather read lengths, completing the first step toward full SBC support.
