# Roadmap: SAM/BAM/CRAM Support in `reads`

## Phase 1: Support Unaligned SAM/BAM/CRAM in rasusa (Done)
- [x] 01-01-PLAN.md — Prototype SBC reader and implement read_lengths.
- [x] 01-02-PLAN.md — Implement filter_reads_into and writing for SBC.
- [x] 01-03-PLAN.md — Paired-end support via SEGMENTED flag and CLI update.

## Phase 2: Refactoring & Abstraction (Done)
- [x] Refactor `src/reads.rs` to use a more generic `RecordSource` abstraction instead of hardcoded `Fastx`.
- [x] Implement `RecordSource` for `Fastx` (current behavior).
- [x] Implement `RecordSource` for `Alignment` (SAM/BAM/CRAM via `noodles`).

## Phase 3: Implementation (Done)
- [x] Implement `read_lengths` for SAM/BAM/CRAM.
- [x] Implement `filter_reads_into` for SAM/BAM/CRAM.
- [x] Support paired-end logic for SAM/BAM/CRAM.
- [x] Update CLI to accept SAM/BAM/CRAM extensions in the `reads` subcommand.

## Phase 4: Validation & Testing (Done)
- [x] Add unit tests for SAM/BAM/CRAM record length gathering.
- [x] Add integration tests using small SAM/BAM/CRAM test cases.
- [x] Verify paired-end integrity for single-file BAM inputs.
- [x] Ensure output format matches input or `-O` flag.

## Phase 5: Documentation & Cleanup
- [ ] Update README with new format support.
- [ ] Final code review and performance check.
