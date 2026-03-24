# Technology Stack

**Analysis Date:** 2025-01-24

## Languages

**Primary:**
- Rust (Edition 2018) - Core application logic and CLI implementation.

**Secondary:**
- Bash - Used in `install/install.sh` and CI scripts.
- Dockerfile - Used for containerization.

## Runtime

**Environment:**
- Native (compiled binary)

**Package Manager:**
- Cargo (Rust's package manager)
- Lockfile: `Cargo.lock` present.

## Frameworks

**Core:**
- `clap` (v4.5.60) - Command-line argument parsing with derive macros.
- `noodles` (v0.104.0) - High-level bioinformatics library for SAM, BAM, and CRAM formats.
- `needletail` (v0.5.1) - Fast sequence file parsing for FASTQ/FASTA formats.

**Testing:**
- `assert_cmd` (v2.0.17) - Testing CLI commands and their output.
- `predicates` (v3.1.4) - Boolean predicates for assertion logic.
- `tempfile` (v3.25.0) - Managing temporary files during tests.
- `cargo-tarpaulin` - Code coverage tool for Rust.

**Build/Dev:**
- `just` - Command runner used for task automation (`justfile`).
- `cross` - Used for cross-compilation to multiple targets in CI.
- `cocogitto` (cog) - Used for conventional commits and versioning (`cog.toml`).

## Key Dependencies

**Critical:**
- `noodles` & `noodles-util` - Essential for alignment file (BAM/SAM) processing.
- `needletail` - Essential for sequence file (FASTQ) processing.
- `rand` & `rand_pcg` - Used for random subsampling logic.
- `niffler` - Provides transparent decompression for various compression formats.

**Infrastructure:**
- `anyhow` & `thiserror` - Standardized error handling patterns.
- `log` & `env_logger` - Logging infrastructure.

## Configuration

**Environment:**
- Configured via command-line arguments (parsed by `clap`).
- `RUST_LOG` environment variable for controlling log verbosity.

**Build:**
- `Cargo.toml` - Main project configuration.
- `justfile` - Task definitions for testing, linting, and building.
- `.release-please-config.json` - Configuration for automated releases.

## Platform Requirements

**Development:**
- Rust toolchain (stable).
- `just` (optional, but recommended).
- `docker` (for container builds).

**Production:**
- Compiled binaries available for:
    - `x86_64-unknown-linux-musl`
    - `x86_64-unknown-linux-gnu`
    - `aarch64-unknown-linux-gnu`
    - `aarch64-unknown-linux-musl`
    - `x86_64-apple-darwin`
    - `aarch64-apple-darwin`
- Docker image available on `ghcr.io`.

---

*Stack analysis: 2025-01-24*
