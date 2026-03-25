# Phase 2: Refactoring & Abstraction - Research

**Researched:** 2026-03-25
**Domain:** Rust Traits, Dynamic Dispatch, Refactoring
**Confidence:** HIGH

## Summary

The current implementation in `src/fastx.rs` relies on branching logic (`if matches!(ext.as_str(), "sam" | "bam" | "cram")`) within `Fastx::read_lengths` and `Fastx::filter_reads_into`. This violates the Single Responsibility Principle and Open/Closed Principle. 

To resolve this, we will introduce a `RecordSource` trait. We will extract the SAM/BAM/CRAM logic into a new `AlignmentSource` struct, and keep the FASTA/FASTQ logic in the `Fastx` struct (or rename it to `FastxSource`). In `src/reads.rs`, we will use dynamic dispatch (`Box<dyn RecordSource>`) to abstract away the file format differences, choosing the correct implementation at runtime via a factory function.

**Primary recommendation:** Introduce a `RecordSource` trait with object-safe methods and a factory function to instantiate `Box<dyn RecordSource>` based on the file extension.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Rust standard library (`std::boxed::Box`) | N/A | Dynamic Dispatch | Allows `reads.rs` to process abstract sources without generic bloat. |
| `std::io::Write` | N/A | Object-safe Writer | Transitioning generic `<T: Write>` to `&mut dyn Write` ensures trait object-safety. |

## Architecture Patterns

### Recommended Project Structure
```text
src/
├── reads.rs       # Uses Box<dyn RecordSource> for generic read processing
├── fastx.rs       # (Or source.rs) Defines the RecordSource trait, FastxSource, and AlignmentSource
├── alignment.rs   # Existing alignment subsampling logic
```

### Pattern 1: Strategy / Polymorphism (Object-Safe Trait)
**What:** Define a trait for format-agnostic sequence reading.
**When to use:** When processing multiple file formats that share a common abstract behavior (e.g. read lengths, filtering).
**Example:**
```rust
pub trait RecordSource {
    fn read_lengths(&self) -> Result<Vec<u32>, FastxError>;
    
    // Note: Use `&mut dyn Write` instead of a generic parameter `<T: Write>` 
    // to keep the trait object-safe, allowing us to create a `Box<dyn RecordSource>`.
    fn filter_reads_into(
        &self,
        reads_to_keep: &[bool],
        nb_reads_keep: usize,
        write_to: &mut dyn std::io::Write,
        output_format: Option<noodles_util::alignment::io::Format>,
    ) -> Result<usize, FastxError>;
}
```

### Pattern 2: Factory Function
**What:** A factory function to determine the format and return the appropriate struct wrapped in a `Box`.
**When to use:** At the entry point in `reads.rs` to encapsulate format sniffing.
**Example:**
```rust
pub fn determine_record_source(path: &Path) -> Box<dyn RecordSource> {
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_lowercase();
    if matches!(ext.as_str(), "sam" | "bam" | "cram") {
        Box::new(AlignmentSource::from_path(path))
    } else {
        Box::new(Fastx::from_path(path))
    }
}
```

### Anti-Patterns to Avoid
- **Generic Trait Methods (`fn filter_reads_into<T: Write>`)**: This makes the trait non-object-safe, preventing the use of `Box<dyn RecordSource>`. Always use trait objects (`&mut dyn std::io::Write`) for parameters in a polymorphic trait.
- **Putting formatting / writer creation into the Source trait**: `Fastx::create` (which creates the output handle) is not part of "reading" sequences. It should be extracted into a standalone utility function (e.g. `create_output_writer(path, compression_level, format)`) rather than being a method on `RecordSource`.

## Runtime State Inventory

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| Stored data | None | None |
| Live service config | None | None |
| OS-registered state | None | None |
| Secrets/env vars | None | None |
| Build artifacts | Rust targets (`target/`) | Clean/rebuild after refactoring |

## Common Pitfalls

### Pitfall 1: Object-Safety Violations
**What goes wrong:** The Rust compiler emits an error like "the trait `RecordSource` cannot be made into an object".
**Why it happens:** One of the methods in the trait has generic type parameters (e.g. `<T: Write>`) or returns `Self`.
**How to avoid:** Use `&mut dyn std::io::Write` instead of `<T: Write>` for the output writer parameter in `filter_reads_into`.

### Pitfall 2: Output File Handle Creation
**What goes wrong:** `reads.rs` currently relies on `out_fastx.create(...)` to build the `output_handle`.
**Why it happens:** Because `Fastx` currently handles both input and output file logic.
**How to avoid:** Decouple output creation from the `Fastx` struct. Create a freestanding `pub fn create_output_writer(...) -> Result<Box<dyn std::io::Write>, FastxError>` function to replace `Fastx::create`.

## Code Examples

### 1. Refactoring `reads.rs` entry points
```rust
// Old:
let input_fastx = Fastx::from_path(&self.input[0]);

// New:
let input_source = determine_record_source(&self.input[0]);
// input_source is a Box<dyn RecordSource>
let mut read_lengths = input_source.read_lengths().context(...)?;
```

### 2. Standalone Output Writer Creation
```rust
pub fn create_output_writer(
    path: &Path,
    compression_lvl: Option<niffler::compression::Level>,
    compression_fmt: Option<niffler::compression::Format>,
) -> Result<Box<dyn std::io::Write>, FastxError> {
    let file = std::fs::File::create(path).map_err(|source| FastxError::CreateError { source })?;
    let file_handle = Box::new(std::io::BufWriter::new(file));
    let fmt = compression_fmt.unwrap_or_else(|| niffler::Format::from_path(path));
    
    let compression_lvl = compression_lvl.unwrap_or(match fmt {
        niffler::compression::Format::Gzip => niffler::compression::Level::Six,
        niffler::compression::Format::Bzip => niffler::compression::Level::Nine,
        niffler::compression::Format::Lzma => niffler::compression::Level::Six,
        niffler::compression::Format::Zstd => niffler::compression::Level::Three,
        _ => niffler::compression::Level::Zero,
    });
    niffler::get_writer(file_handle, fmt, compression_lvl)
        .map_err(FastxError::CompressOutputError)
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Branching on extension in `Fastx` methods | `RecordSource` trait + Factory | Phase 2 | Clean separation of concerns, easier to test and extend. |

## Open Questions

1. **Where should `RecordSource` live?**
   - Recommendation: Create a new module `src/source.rs` containing `RecordSource`, `AlignmentSource`, and `FastxSource`, leaving `src/fastx.rs` for shared logic, OR rename `src/fastx.rs` to something more generic if appropriate. Given typical refactoring scopes, putting the trait in `src/fastx.rs` (which already defines `FastxError`) and just splitting the structs inside the same file might be simplest.

## Sources

### Primary (HIGH confidence)
- Codebase context (`src/fastx.rs`, `src/reads.rs`)

### Secondary (MEDIUM confidence)
- Rust Reference: [Trait Objects](https://doc.rust-lang.org/reference/items/traits.html#object-safety)

## Metadata

**Confidence breakdown:**
- Architecture: HIGH - Standard Rust dynamic dispatch refactoring.
- Pitfalls: HIGH - Generic traits preventing object safety is a very common Rust pitfall.

**Research date:** 2026-03-25
**Valid until:** Indefinite (internal codebase refactoring)