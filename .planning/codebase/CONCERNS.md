# Codebase Concerns

**Analysis Date:** 2024-12-19

## Tech Debt

**Large File and Complex Logic:**
- Issue: `src/alignment.rs` is excessively large (nearly 1900 lines) and handles multiple complex subsampling strategies (Stream and Fetch) for different alignment formats (SAM/BAM/CRAM).
- Files: `src/alignment.rs`
- Impact: High maintenance cost, difficult to test individual components, and increased risk of regressions when modifying logic.
- Fix approach: Refactor `src/alignment.rs` into a module with separate files for each strategy (e.g., `src/alignment/stream.rs`, `src/alignment/fetch.rs`) and shared utilities.

**In-Memory Read Length Storage:**
- Issue: For FASTA/FASTQ subsampling, the tool reads the entire file twice. The first pass stores all read lengths in a `Vec<u32>`.
- Files: `src/fastx.rs`, `src/reads.rs`
- Impact: High memory consumption for large datasets (e.g., billions of reads). Significant performance hit due to double file I/O.
- Fix approach: Implement a single-pass subsampling algorithm (e.g., Reservoir Sampling) or use a more memory-efficient way to track reads if possible.

**Direct Process Exit in Library-like Code:**
- Issue: The `run` method in `src/reads.rs` calls `std::process::exit(1)` directly when paired-end files have mismatched read counts.
- Files: `src/reads.rs`
- Impact: Makes the code harder to test and prevents graceful handling of errors if this logic were ever used as a library.
- Fix approach: Return a `Result` with a custom error variant instead of exiting the process.

## Performance Bottlenecks

**Double Pass on FASTX Files:**
- Problem: `reads` subcommand reads the input file once to get lengths and once to filter.
- Files: `src/reads.rs`, `src/fastx.rs`
- Cause: Current subsampling logic requires knowing all lengths (or total count) to generate indices before filtering.
- Improvement path: Switch to a reservoir sampling or another single-pass approach that doesn't require prior knowledge of the total number of reads.

**Owned Record Conversion:**
- Problem: Alignment records are frequently converted from borrowed `Record` to owned `RecordBuf`.
- Files: `src/alignment.rs`
- Cause: Simplifies logic for storage in heaps and caches, but at the cost of memory and CPU cycles.
- Improvement path: Use borrowed records where possible, or use more efficient storage mechanisms.

## Fragile Areas

**Paired-End Survivor Tracking:**
- Files: `src/alignment.rs`
- Why fragile: Uses `HashSet<Vec<u8>>` to store all survivor read names to find their mates.
- Safe modification: Be careful with memory usage here; millions of read names can easily consume gigabytes of RAM.
- Test coverage: Gaps in testing very large paired-end alignment files where this set would grow significantly.

## Scaling Limits

**32-bit Index Truncation:**
- Current capacity: $2^{32} - 1$ reads (approx. 4.2 billion).
- Limit: `SubSampler::shuffled_indices` casts `usize` to `u32`.
- Files: `src/subsampler.rs`
- Scaling path: Change indices to `u64` or `usize` to support larger datasets, and update the subsampling logic to avoid full index vector generation.

**Memory-Bound FASTX Subsampling:**
- Current capacity: Limited by available RAM to store `Vec<u32>` of lengths.
- Limit: Approximately 4 bytes per read. 1 billion reads = 4GB RAM just for lengths.
- Scaling path: Single-pass sampling or disk-backed index.

## Test Coverage Gaps

**Large Dataset Edge Cases:**
- What's not tested: Performance and memory behavior with multi-gigabyte files or files with >4B reads.
- Files: `src/fastx.rs`, `src/alignment.rs`, `src/subsampler.rs`
- Risk: Tool may crash or produce incorrect results (due to truncation) on very large genomic datasets.
- Priority: Medium

---

*Concerns audit: 2024-12-19*
