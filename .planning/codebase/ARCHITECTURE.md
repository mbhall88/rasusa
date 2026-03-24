# Architecture

**Analysis Date:** 2025-02-11

## Pattern Overview

**Overall:** Command-line tool with subcommand-based execution using a Trait-based polymorphic approach.

**Key Characteristics:**
- **Trait-based Dispatch:** Uses a `Runner` trait to abstract subcommand execution logic.
- **Library-backed I/O:** Leverages specialized bioinformatics libraries (`needletail` for reads, `noodles` for alignments) for robust format support.
- **Deterministic Randomness:** Uses PCG-based random number generation (`rand_pcg`) to allow reproducible subsampling via seeds.

## Layers

**CLI Layer:**
- Purpose: Handles argument parsing, validation, and user interface.
- Location: `src/cli.rs`
- Contains: `clap` structs, custom parsers for biological units (GenomeSize, Coverage), and error definitions.
- Depends on: `clap`, `regex`.
- Used by: `src/main.rs`.

**Command Layer (Subcommands):**
- Purpose: Orchestrates the subsampling process for different data types.
- Location: `src/reads.rs`, `src/alignment.rs`
- Contains: Implementations of the `Runner` trait for `Reads` and `Alignment` subcommands.
- Depends on: CLI Layer, I/O Wrappers, Subsampling Core.
- Used by: `src/main.rs`.

**I/O Wrapper Layer:**
- Purpose: Provides a simplified interface for reading/writing biological file formats with transparent compression.
- Location: `src/fastx.rs`
- Contains: `Fastx` struct for FASTA/FASTQ handling.
- Depends on: `needletail`, `niffler`.
- Used by: `src/reads.rs`.

**Subsampling Core:**
- Purpose: Implements the mathematical and algorithmic logic for selecting a subset of data.
- Location: `src/subsampler.rs`, `src/alignment.rs`
- Contains: `SubSampler` (for reads) and `stream`/`fetching` methods (for alignments).
- Depends on: `rand`, `rand_pcg`.
- Used by: Command Layer.

## Data Flow

**Read Subsampling Flow:**

1. `src/main.rs` parses arguments into a `Reads` struct.
2. `Reads::run()` validates inputs and outputs.
3. `Fastx` is used to pre-scan input files for read lengths.
4. `SubSampler` generates a list of indices to keep based on the target (coverage, bases, or count).
5. `Fastx` performs a second pass, streaming selected reads to the output.

**Alignment Subsampling Flow:**

1. `src/main.rs` parses arguments into an `Alignment` struct.
2. `Alignment::run()` selects between `Stream` or `Fetch` strategy.
3. **Stream Strategy:** Performs a single-pass sweep-line algorithm using a Max-Heap to maintain target coverage at every position.
4. **Fetch Strategy:** Uses index-based random access to sample regions across the genome.
5. For paired-end data, a second pass (`recover_mates`) identifies mate pairs for selected reads.

**State Management:**
- Stateless between runs.
- In-memory state is minimized using streaming patterns (sweep-line for alignments, index-based selection for reads).

## Key Abstractions

**Runner Trait:**
- Purpose: Unified interface for executing subcommands.
- Examples: `src/main.rs`
- Pattern: Strategy Pattern.

**Fastx Struct:**
- Purpose: Unified handling of FASTA/FASTQ files with various compression formats.
- Examples: `src/fastx.rs`
- Pattern: Wrapper/Facade.

**ScoredRead Struct:**
- Purpose: Internal representation of an alignment record with a random priority key for the sweep-line algorithm.
- Examples: `src/alignment.rs`

## Entry Points

**main function:**
- Location: `src/main.rs`
- Triggers: User execution via CLI.
- Responsibilities: CLI parsing, Logger initialization, Subcommand dispatch via `Box<dyn Runner>`.

## Error Handling

**Strategy:** Result-based error propagation using `anyhow` for top-level errors and `thiserror` for library-specific error types.

**Patterns:**
- `CliError`: Custom enum in `src/cli.rs` for validation errors.
- `FastxError`: Custom enum in `src/fastx.rs` for I/O and parsing errors.
- Top-level `anyhow::Result` in `main.rs` for clean user-facing error messages.

## Cross-Cutting Concerns

**Logging:** Uses `env_logger` and `log` crate. Verbosity controlled via `-v` flag.
**Validation:** Complex validation logic for biological units (e.g., parsing `4.3kb` to bytes) in `src/cli.rs`.
**Authentication:** Not applicable (local file processing).

---

*Architecture analysis: 2025-02-11*
