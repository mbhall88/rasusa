//! Single source of truth for file-format/compression detection and the mappings derived from
//! it. Before this module existed, extension->format detection, the compression-level defaults,
//! and the `OutputFormat`->noodles `Format` mapping were each duplicated across `reads.rs`,
//! `fastx.rs`, and `alignment/io.rs`, with no guarantee the copies agreed.

use std::path::Path;

use niffler::compression;
use noodles_util::alignment::io::Format;

use crate::cli::OutputFormat;

/// What [`crate::source::RecordSource::filter_reads_into`] should encode records as. Threading
/// this through (rather than a raw `Option<Format>` + `bool` pair) means implementations that
/// can never receive alignment output (e.g. [`crate::fastx::Fastx`]) don't need to reason about
/// the noodles `Format` type at all.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputEncoding {
    /// Write as SAM/BAM/CRAM.
    Alignment(Format),
    /// Write as FASTA (`fasta: true`) or FASTQ (`fasta: false`).
    Fastx { fasta: bool },
}

/// Infers an alignment [`Format`] from a path's extension (`.sam`/`.bam`/`.cram`).
pub fn infer_format_from_path(path: &Path) -> Option<Format> {
    match path.extension().and_then(|ext| ext.to_str()) {
        Some("sam") => Some(Format::Sam),
        Some("bam") => Some(Format::Bam),
        Some("cram") => Some(Format::Cram),
        _ => None,
    }
}

/// Infers an alignment [`Format`] from a single character/word (case-insensitive), for CLI
/// value parsing.
pub fn infer_format_from_char(c: &str) -> Result<Format, String> {
    match c.to_ascii_lowercase().as_str() {
        "s" | "sam" => Ok(Format::Sam),
        "b" | "bam" => Ok(Format::Bam),
        "c" | "cram" => Ok(Format::Cram),
        _ => Err(String::from("Invalid output format. Please use 's', 'b', or 'c' to specify SAM, BAM, or CRAM format, respectively")),
    }
}

/// Whether `path` names a FASTA file, based on its extension (a trailing compression extension
/// such as `.gz` is stripped first).
pub fn is_fasta_extension(path: &Path) -> bool {
    let mut p = path.to_path_buf();
    if let Some(ext) = p.extension().map(|e| e.to_string_lossy().to_lowercase()) {
        if ["gz", "bz2", "xz", "zst", "bz"].contains(&ext.as_str()) {
            p = p.with_extension("");
        }
    }
    if let Some(ext) = p.extension().map(|e| e.to_string_lossy().to_lowercase()) {
        return ext == "fasta" || ext == "fa";
    }
    false
}

/// Whether output should be written as FASTA, given an explicit `--output-format` (if any) and a
/// path to fall back to extension-based detection with when no explicit choice was made.
pub fn is_fasta_output(output_format: Option<&OutputFormat>, fallback_path: &Path) -> bool {
    match output_format {
        Some(OutputFormat::Fasta) => true,
        Some(OutputFormat::Fastq) => false,
        _ => is_fasta_extension(fallback_path),
    }
}

/// Maps a CLI-level `--output-format` selection to the noodles alignment [`Format`] it
/// corresponds to. Returns `None` for FASTA/FASTQ selections or when no `--output-format` was
/// given - callers fall back to extension-based inference ([`infer_format_from_path`]) in that
/// case.
pub fn output_alignment_format(output_format: Option<&OutputFormat>) -> Option<Format> {
    match output_format {
        Some(OutputFormat::Sam) => Some(Format::Sam),
        Some(OutputFormat::Bam) => Some(Format::Bam),
        Some(OutputFormat::Cram) => Some(Format::Cram),
        Some(OutputFormat::Fasta) | Some(OutputFormat::Fastq) | None => None,
    }
}

/// Infers a compression format from a path's extension (e.g. `.gz` -> Gzip), falling back to
/// [`compression::Format::No`] for an unrecognised or absent extension.
pub fn infer_compression_format(path: &Path) -> compression::Format {
    match path.extension().and_then(|ext| ext.to_str()) {
        Some("gz") => compression::Format::Gzip,
        Some("bz") | Some("bz2") => compression::Format::Bzip,
        Some("lzma") => compression::Format::Lzma,
        _ => compression::Format::No,
    }
}

/// The default compression level rasusa uses for a given compression format, when the user
/// hasn't requested a specific level.
pub fn default_compression_level(fmt: compression::Format) -> compression::Level {
    match fmt {
        compression::Format::Gzip => compression::Level::Six,
        compression::Format::Bzip => compression::Level::Nine,
        compression::Format::Lzma => compression::Level::Six,
        compression::Format::Zstd => compression::Level::Three,
        _ => compression::Level::Zero,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_infer_compression_format() {
        assert_eq!(
            infer_compression_format(Path::new("foo.gz")),
            compression::Format::Gzip
        );
        assert_eq!(
            infer_compression_format(Path::new("baz")),
            compression::Format::No
        );
        assert_eq!(
            infer_compression_format(Path::new("baz.fq")),
            compression::Format::No
        );
        assert_eq!(
            infer_compression_format(Path::new("baz.fq.bz2")),
            compression::Format::Bzip
        );
        assert_eq!(
            infer_compression_format(Path::new("baz.fq.bz")),
            compression::Format::Bzip
        );
        assert_eq!(
            infer_compression_format(Path::new("baz.fq.lzma")),
            compression::Format::Lzma
        );
    }

    #[test]
    fn test_infer_format() {
        let fmt = infer_format_from_path(Path::new("file.sam")).unwrap();
        assert!(matches!(fmt, Format::Sam));

        let fmt = infer_format_from_path(Path::new("file.bam")).unwrap();
        assert!(matches!(fmt, Format::Bam));

        let fmt = infer_format_from_path(Path::new("file.cram")).unwrap();
        assert!(matches!(fmt, Format::Cram));

        let fmt = infer_format_from_path(Path::new("file.txt"));
        assert!(fmt.is_none());
    }

    #[test]
    fn test_infer_format_from_char() {
        let fmt = infer_format_from_char("s").unwrap();
        assert!(matches!(fmt, Format::Sam));

        let fmt = infer_format_from_char("b").unwrap();
        assert!(matches!(fmt, Format::Bam));

        let fmt = infer_format_from_char("c").unwrap();
        assert!(matches!(fmt, Format::Cram));

        let fmt = infer_format_from_char("x");
        assert!(fmt.is_err());

        let fmt = infer_format_from_char("sam").unwrap();
        assert!(matches!(fmt, Format::Sam));

        let fmt = infer_format_from_char("B").unwrap();
        assert!(matches!(fmt, Format::Bam));

        let fmt = infer_format_from_char("CRAM").unwrap();
        assert!(matches!(fmt, Format::Cram));

        let fmt = infer_format_from_char("x");
        assert!(fmt.is_err());
    }

    #[test]
    fn test_is_fasta_extension() {
        assert!(is_fasta_extension(Path::new("reads.fasta")));
        assert!(is_fasta_extension(Path::new("reads.fa")));
        assert!(is_fasta_extension(Path::new("reads.fa.gz")));
        assert!(!is_fasta_extension(Path::new("reads.fastq")));
        assert!(!is_fasta_extension(Path::new("reads.fq.gz")));
    }

    #[test]
    fn test_is_fasta_output_explicit_overrides_fallback_path() {
        assert!(is_fasta_output(
            Some(&OutputFormat::Fasta),
            Path::new("reads.fastq")
        ));
        assert!(!is_fasta_output(
            Some(&OutputFormat::Fastq),
            Path::new("reads.fasta")
        ));
    }

    #[test]
    fn test_is_fasta_output_falls_back_to_path_extension() {
        assert!(is_fasta_output(None, Path::new("reads.fasta")));
        assert!(!is_fasta_output(None, Path::new("reads.fastq")));
        // Sam/Bam/Cram aren't Fasta/Fastq, so this also falls back to the path extension - the
        // resulting bool is only ever consulted by callers when the *output encoding* is Fastx.
        assert!(is_fasta_output(
            Some(&OutputFormat::Sam),
            Path::new("reads.fasta")
        ));
    }

    #[test]
    fn test_output_alignment_format() {
        assert_eq!(
            output_alignment_format(Some(&OutputFormat::Sam)),
            Some(Format::Sam)
        );
        assert_eq!(
            output_alignment_format(Some(&OutputFormat::Bam)),
            Some(Format::Bam)
        );
        assert_eq!(
            output_alignment_format(Some(&OutputFormat::Cram)),
            Some(Format::Cram)
        );
        assert_eq!(output_alignment_format(Some(&OutputFormat::Fasta)), None);
        assert_eq!(output_alignment_format(Some(&OutputFormat::Fastq)), None);
        assert_eq!(output_alignment_format(None), None);
    }

    #[test]
    fn test_default_compression_level() {
        assert_eq!(
            default_compression_level(compression::Format::Gzip),
            compression::Level::Six
        );
        assert_eq!(
            default_compression_level(compression::Format::Bzip),
            compression::Level::Nine
        );
        assert_eq!(
            default_compression_level(compression::Format::Lzma),
            compression::Level::Six
        );
        assert_eq!(
            default_compression_level(compression::Format::Zstd),
            compression::Level::Three
        );
        assert_eq!(
            default_compression_level(compression::Format::No),
            compression::Level::Zero
        );
    }
}
