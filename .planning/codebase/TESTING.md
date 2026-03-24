# Testing Patterns

**Analysis Date:** 2025-03-24

## Test Framework

**Runner:**
- Rust's built-in `cargo test`.
- Orchestrated by `just` in CI: `just test` (`cargo test -v --all-targets --no-fail-fast`).

**Assertion Library:**
- Standard library `assert!`, `assert_eq!`.
- `predicates` for complex assertions on command output.
- `assert_cmd` for binary testing.

**Run Commands:**
```bash
cargo test              # Run all tests
just test               # Run all tests (via just)
cargo test --doc        # Run documentation tests
just coverage           # Run coverage via tarpaulin
```

## Test File Organization

**Location:**
- **Unit Tests**: Co-located within source files in `src/*.rs`.
- **Integration Tests**: In the `tests/` directory.
- **Doc-tests**: Inline within documentation comments in `src/*.rs`.

**Naming:**
- Unit test modules: `#[cfg(test)] mod tests`.
- Integration test files: `main.rs`, `reproducibility.rs`.

**Structure:**
```
src/
  cli.rs (contains mod tests)
  ...
tests/
  main.rs          # CLI/integration tests
  reproducibility.rs # Specific seed-based reproducibility checks
  cases/           # Test input files (FastQ, BAM, FASTA, etc.)
```

## Test Structure

**Suite Organization:**
```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_name() {
        // setup
        // assert
    }
}
```

**Patterns:**
- **Unit Testing**: Tests individual functions/types in the same file where they are defined.
- **CLI Testing**: Uses `assert_cmd::Command::cargo_bin(BIN)?` to execute the binary and check behavior.
- **Reproducibility**: Checks output for deterministic subsampling by comparing against known baseline files for a given seed.

## Mocking

**Framework:** None.
- Not extensively used as the logic is data-processing centric.

**Patterns:**
- Testing with real data files in `tests/cases/`.

**What to Mock:**
- Not applicable.

**What NOT to Mock:**
- File I/O (often tested with real files or `tempfile`).
- CLI behavior (tested via `assert_cmd`).

## Fixtures and Factories

**Test Data:**
```rust
// In tests/reproducibility.rs
let cases = vec![
    (
        1,
        vec!["@read1", "@read2", ...],
    ),
    ...
];
```

**Location:**
- Test input files live in `tests/cases/`.
- Hardcoded expected outputs for small cases are within the test files.

## Coverage

**Requirements:** `codecov.yml` is present, suggesting coverage is tracked in CI.

**View Coverage:**
```bash
just coverage # uses cargo tarpaulin
```

## Test Types

**Unit Tests:**
- Test logic in `src/cli.rs`, `src/alignment.rs`, etc. Focus on parsing and calculations.

**Integration Tests:**
- In `tests/main.rs`, tests end-to-end execution of the `rasusa` binary with various flags and inputs.

**E2E Tests:**
- Included in integration tests. Checks output file existence and content.

## Common Patterns

**Async Testing:** Not used.

**Error Testing:**
```rust
#[test]
fn input_file_doesnt_exist() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin(BIN)?;
    cmd.args(vec![READS, "file/doesnt/exist.fa", "-g", "5mb", "-c", "20"]);
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("does not exist"));

    Ok(())
}
```

---

*Testing analysis: 2025-03-24*
