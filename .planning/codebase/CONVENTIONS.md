# Coding Conventions

**Analysis Date:** 2025-03-24

## Naming Patterns

**Files:**
- Snake case: `alignment.rs`, `cli.rs`, `fastx.rs`.

**Functions:**
- Snake case: `parse_fraction`, `check_path_exists`.

**Variables:**
- Snake case: `log_lvl`, `log_builder`.

**Types/Traits/Enums:**
- PascalCase: `Cli`, `Commands`, `CliError`, `Runner`, `GenomeSize`, `Coverage`.

**Constants:**
- SCREAMING_SNAKE_CASE: `CITATION`, `BIN`, `READS`.

## Code Style

**Formatting:**
- Standard Rust formatting using `rustfmt`.
- Enforced in CI via `cargo fmt --all -- --check`.

**Linting:**
- `clippy` is used for linting.
- Configuration in `.clippy.toml`: `cognitive-complexity-threshold = 30`.
- Enforced in CI via `just lint` (`cargo clippy --all-features --all-targets -- -D warnings`).

## Import Organization

**Order:**
1. Standard library imports (`std::*`)
2. External crate imports (`clap::*`, `log::*`)
3. Internal module imports (`crate::*`)

**Path Aliases:**
- Not extensively used. `pub use crate::cli::Cli;` is used for re-exporting in `src/main.rs`.

## Error Handling

**Patterns:**
- **Library/CLI Errors**: Uses `thiserror` crate to define custom error enums with clear messages. Example: `CliError` in `src/cli.rs`.
- **Top-level/Application Errors**: Uses `anyhow::Result` for application-level result handling in `main` and high-level traits like `Runner`.
- **Validation**: Helper functions like `check_path_exists` return `Result<PathBuf, String>`.

## Logging

**Framework:** `log` crate with `env_logger` as the backend.

**Patterns:**
- Verbosity controlled via `-v/--verbose` CLI flag in `src/main.rs`.
- `debug!`, `info!`, etc. used throughout the codebase for trace and status information.

## Comments

**When to Comment:**
- Complexity: Explain non-obvious logic (e.g., regex in `src/cli.rs`).
- Public API: All public types and functions generally have documentation comments.

**JSDoc/TSDoc:**
- Triple-slash `///` comments used for documentation.
- Triple-slash comments often include `# Example` sections that serve as doc-tests.

## Function Design

**Size:**
- Generally small and focused. `clippy` enforces cognitive complexity limits.

**Parameters:**
- Uses `AsRef<OsStr>` or `AsRef<Path>` for file path arguments to be flexible.

**Return Values:**
- Returns `Result` for operations that can fail.
- Returns `self` in some trait implementations for chaining or conversion.

## Module Design

**Exports:**
- `src/main.rs` re-exports key types: `pub use crate::cli::Cli;`, `pub use crate::fastx::Fastx;`.

**Barrel Files:**
- `src/main.rs` acts as a barrel file for internal modules.

---

*Convention analysis: 2025-03-24*
