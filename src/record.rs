use crate::fastx::FastxError;
use std::io::Write;

/// Writes a single FASTA or FASTQ record to `write_to`, using `line_ending` (e.g. `b"\n"` or
/// `b"\r\n"`) to terminate every line.
///
/// Writes FASTQ (`@name<eol>seq<eol>+<eol>qual<eol>`) when `qual` is present and `fasta` is
/// `false`; otherwise writes FASTA (`>name<eol>seq<eol>`). This mirrors the two existing record
/// sources: a FASTA/Q record carries no quality scores, and an alignment record with no quality
/// scores is represented the same way.
pub fn write_fastx_record(
    write_to: &mut dyn Write,
    name: &[u8],
    seq: &[u8],
    qual: Option<&[u8]>,
    fasta: bool,
    line_ending: &[u8],
) -> Result<(), FastxError> {
    let map_err = |source: std::io::Error| FastxError::WriteError {
        source: anyhow::Error::from(source),
    };

    match qual {
        Some(qual) if !fasta => {
            write_to.write_all(b"@").map_err(map_err)?;
            write_to.write_all(name).map_err(map_err)?;
            write_to.write_all(line_ending).map_err(map_err)?;
            write_to.write_all(seq).map_err(map_err)?;
            write_to.write_all(line_ending).map_err(map_err)?;
            write_to.write_all(b"+").map_err(map_err)?;
            write_to.write_all(line_ending).map_err(map_err)?;
            write_to.write_all(qual).map_err(map_err)?;
            write_to.write_all(line_ending).map_err(map_err)?;
        }
        _ => {
            write_to.write_all(b">").map_err(map_err)?;
            write_to.write_all(name).map_err(map_err)?;
            write_to.write_all(line_ending).map_err(map_err)?;
            write_to.write_all(seq).map_err(map_err)?;
            write_to.write_all(line_ending).map_err(map_err)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn writes_fasta_when_fasta_requested_even_with_quality() {
        let mut out: Vec<u8> = Vec::new();
        write_fastx_record(&mut out, b"read1", b"ACGT", Some(b"!!!!"), true, b"\n").unwrap();
        assert_eq!(out, b">read1\nACGT\n");
    }

    #[test]
    fn writes_fasta_when_no_quality_present() {
        let mut out: Vec<u8> = Vec::new();
        write_fastx_record(&mut out, b"read1", b"ACGT", None, false, b"\n").unwrap();
        assert_eq!(out, b">read1\nACGT\n");
    }

    #[test]
    fn writes_fastq_when_quality_present_and_fastq_requested() {
        let mut out: Vec<u8> = Vec::new();
        write_fastx_record(&mut out, b"read1", b"ACGT", Some(b"!!!!"), false, b"\n").unwrap();
        assert_eq!(out, b"@read1\nACGT\n+\n!!!!\n");
    }

    #[test]
    fn writes_fastq_with_windows_line_endings() {
        let mut out: Vec<u8> = Vec::new();
        write_fastx_record(&mut out, b"read1", b"ACGT", Some(b"!!!!"), false, b"\r\n").unwrap();
        assert_eq!(out, b"@read1\r\nACGT\r\n+\r\n!!!!\r\n");
    }
}
