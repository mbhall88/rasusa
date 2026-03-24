# Codebase Structure

**Analysis Date:** 2025-02-11

## Directory Layout

```
rasusa/
├── src/
│   ├── main.rs         # Application entry point and subcommand dispatch
│   ├── cli.rs          # CLI argument parsing and custom unit parsers
│   ├── reads.rs        # Read subsampling subcommand implementation
│   ├── alignment.rs    # Alignment subsampling subcommand implementation
│   ├── fastx.rs        # FASTA/FASTQ I/O wrapper
│   ├── subsampler.rs   # Core read subsampling algorithm
│   └── reads.rs        # Collection of read sets
├── tests/
│   ├── cases/          # Test data for integration tests
│   ├── main.rs         # Integration tests
│   └── reproducibility.rs # Tests for consistent seeding
├── .github/
│   └── workflows/      # CI/CD pipelines
├── paper/              # Source for JOSS publication
└── img/                # Project assets (logos, etc.)
```

## Directory Purposes

**src/:**
- Purpose: Primary application source code.
- Contains: Rust source files.
- Key files: `src/main.rs`, `src/cli.rs`.

**tests/:**
- Purpose: Integration and reproducibility testing.
- Contains: Rust test files and a subdirectory for test cases.
- Key files: `tests/main.rs`, `tests/cases/`.

**paper/:**
- Purpose: Scientific publication materials.
- Contains: LaTeX/Markdown sources for the JOSS paper.

## Key File Locations

**Entry Points:**
- `src/main.rs`: Dispatches to `Reads` or `Alignment` subcommands.

**Configuration:**
- `Cargo.toml`: Dependency management and project metadata.
- `.clippy.toml`: Clippy lint configuration.

**Core Logic:**
- `src/subsampler.rs`: Logic for selecting reads to keep.
- `src/alignment.rs`: Logic for alignment subsampling (Stream and Fetch strategies).

**Testing:**
- `tests/main.rs`: Comprehensive CLI integration tests.
- `src/*.rs` (internal): Unit tests co-located within source files.

## Naming Conventions

**Files:**
- Snake case for Rust files: `alignment.rs`, `subsampler.rs`.

**Directories:**
- Snake case or simple names: `tests/cases/`, `paper/`.

## Where to Add New Code

**New Feature (Subcommand):**
1. Define the command struct in `src/cli.rs`.
2. Implement the `Runner` trait in a new file (e.g., `src/new_command.rs`).
3. Add a variant to `Commands` enum in `src/cli.rs`.
4. Update `main` in `src/main.rs` to handle the new variant.

**New Biological Unit/Metric:**
- Add to `src/cli.rs` as a new struct with `FromStr` implementation.

**Utilities:**
- Shared helper functions should go into `src/cli.rs` or a dedicated `utils.rs` (if one is created).

## Special Directories

**tests/cases/:**
- Purpose: Contains compressed and indexed biological data files (BAM, FASTQ, FAIDX).
- Generated: No (checked into version control).
- Committed: Yes.

---

*Structure analysis: 2025-02-11*
