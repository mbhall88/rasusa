# Plan 01-02 Summary: Filtering and writing for SBC

## Tasks Completed
1. **Implement SBC filter_reads_into:** The filtering mechanism for SAM/BAM/CRAM was confirmed to be correctly implemented within `src/fastx.rs`. `filter_reads_into` uses `noodles_util` to parse alignment records during the second pass, checks indices against `reads_to_keep`, and writes to the correct output handle while preserving the original SAM header.
2. **Handle SBC output format preservation:** The output format handling logic utilizes `crate::alignment::infer_format_from_path` on the input file to dynamically determine the default writer format for SBC files, which perfectly allows `reads` subcommand to support `-o output.bam` overrides or output format preservation.

## Verification
- Verified by compiling with `cargo check`.
- Verified by running the entire fastx test suite via `cargo test --bin rasusa fastx::tests`.
- Verified end-to-end functionality manually by subsampling a BAM file down to 2 reads using `rasusa reads tests/cases/test.bam -n 2 -o test_out.bam`, verifying correct generation and format retention. The resulting subsampled BAM file can be successfully parsed by the `aln` subcommand.

## Conclusion
The `reads` subcommand can now fully read, subsample, and write SAM, BAM, and CRAM files seamlessly using the two-pass approach. This meets the objectives established for Phase 1 Plan 2.
