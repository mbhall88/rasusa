# Roadmap: SAM/BAM/CRAM Support in `reads`

## Phase 1: Research & Prototyping
- [ ] Research `noodles` API for reading unaligned SAM/BAM/CRAM records.
- [ ] Determine how to handle paired-end data in a single BAM file (matching by name or assuming adjacency).
- [ ] Prototype a simple unaligned BAM reader in Rust.

## Phase 2: Refactoring & Abstraction
- [ ] Refactor `src/reads.rs` to use a more generic `RecordSource` abstraction instead of hardcoded `Fastx`.
- [ ] Implement `RecordSource` for `Fastx` (current behavior).
- [ ] Implement `RecordSource` for `Alignment` (SAM/BAM/CRAM via `noodles`).

## Phase 3: Implementation
- [ ] Implement `read_lengths` for SAM/BAM/CRAM.
- [ ] Implement `filter_reads_into` for SAM/BAM/CRAM.
- [ ] Support paired-end logic for SAM/BAM/CRAM.
- [ ] Update CLI to accept SAM/BAM/CRAM extensions in the `reads` subcommand.

## Phase 4: Validation & Testing
- [ ] Add unit tests for SAM/BAM/CRAM record length gathering.
- [ ] Add integration tests using small SAM/BAM/CRAM test cases.
- [ ] Verify paired-end integrity for single-file BAM inputs.
- [ ] Ensure output format matches input or `-O` flag.

## Phase 5: Documentation & Cleanup
- [ ] Update README with new format support.
- [ ] Final code review and performance check.
