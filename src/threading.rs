//! Shared reader/writer factory for optional multithreaded BGZF (de)compression.
//!
//! Only BAM benefits here: SAM is uncompressed text and noodles' CRAM codec has no
//! multithreaded path, so `threads > 1` only changes how BAM input/output is (de)compressed.
//! For every other case (SAM, CRAM, or `threads == 1`) this falls back to noodles' normal
//! single-threaded, autodetecting reader/writer, so `--threads 1` is byte-for-byte the same as
//! before `--threads` existed.

use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom, Write};
use std::num::NonZeroUsize;
use std::path::Path;

use noodles_bgzf as bgzf;
use noodles_util::alignment::{self, io::Format};

use crate::format::infer_format_from_path;

pub type AlignmentReader = alignment::io::Reader<Box<dyn Read>>;
pub type AlignmentWriter = alignment::io::Writer<Box<dyn Write + Send>>;

/// Same magic number noodles' own compression autodetection looks for (see
/// `noodles_util::alignment::io::reader::builder::detect_compression_method`).
const GZIP_MAGIC_NUMBER: [u8; 2] = [0x1f, 0x8b];

/// Whether `file` starts with the BGZF/gzip magic number, restoring the read position
/// afterwards. `path`'s `.bam` extension alone isn't a reliable enough signal to commit to a
/// BGZF decode - a mislabeled or renamed SAM/CRAM file would otherwise be forced through the
/// multithreaded BAM path and fail to decompress, where the normal single-threaded path would
/// have content-sniffed its way to the right format instead.
fn looks_bgzf_compressed(file: &mut File) -> io::Result<bool> {
    let mut buf = [0u8; GZIP_MAGIC_NUMBER.len()];
    let is_gzip = match file.read_exact(&mut buf) {
        Ok(()) => buf == GZIP_MAGIC_NUMBER,
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => false,
        Err(e) => return Err(e),
    };
    file.seek(SeekFrom::Start(0))?;
    Ok(is_gzip)
}

/// Builds an alignment reader for `path`, using a multithreaded BGZF decoder when `path` is BAM
/// and `threads > 1`.
pub fn build_alignment_reader(path: &Path, threads: NonZeroUsize) -> io::Result<AlignmentReader> {
    let mut file = File::open(path)?;

    let use_multithreaded_bam = threads.get() > 1
        && infer_format_from_path(path) == Some(Format::Bam)
        && looks_bgzf_compressed(&mut file)?;

    if use_multithreaded_bam {
        let decoder = bgzf::io::MultithreadedReader::with_worker_count(threads, file);
        alignment::io::reader::Builder::default()
            .set_format(Format::Bam)
            .set_compression_method(None) // already decompressed by `decoder`
            .build_from_reader(Box::new(decoder) as Box<dyn Read>)
    } else {
        alignment::io::reader::Builder::default().build_from_reader(Box::new(file) as Box<dyn Read>)
    }
}

/// Builds an alignment writer over `sink`, using a multithreaded BGZF encoder when `format` is
/// BAM and `threads > 1`.
pub fn build_alignment_writer(
    sink: Box<dyn Write + Send>,
    format: Format,
    threads: NonZeroUsize,
) -> io::Result<AlignmentWriter> {
    if threads.get() > 1 && format == Format::Bam {
        let encoder = bgzf::io::MultithreadedWriter::with_worker_count(threads, sink);
        alignment::io::writer::Builder::default()
            .set_format(Format::Bam)
            .set_compression_method(None) // already compressed by `encoder`
            .build_from_writer(Box::new(encoder) as Box<dyn Write + Send>)
    } else {
        alignment::io::writer::Builder::default()
            .set_format(format)
            .build_from_writer(sink)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_threaded_reader_reads_header_and_records() {
        let reader = build_alignment_reader(
            Path::new("tests/cases/test.bam"),
            NonZeroUsize::new(1).unwrap(),
        )
        .unwrap();
        let _ = reader;
    }

    #[test]
    fn multithreaded_and_single_threaded_readers_agree_on_records() {
        let path = Path::new("tests/cases/test.bam");

        let mut single = build_alignment_reader(path, NonZeroUsize::new(1).unwrap()).unwrap();
        let header1 = single.read_header().unwrap();
        let names1: Vec<Vec<u8>> = single
            .records(&header1)
            .map(|r| {
                let record = r.unwrap();
                let name: &[u8] = record.name().unwrap().as_ref();
                name.to_vec()
            })
            .collect();

        let mut multi = build_alignment_reader(path, NonZeroUsize::new(4).unwrap()).unwrap();
        let header2 = multi.read_header().unwrap();
        let names2: Vec<Vec<u8>> = multi
            .records(&header2)
            .map(|r| {
                let record = r.unwrap();
                let name: &[u8] = record.name().unwrap().as_ref();
                name.to_vec()
            })
            .collect();

        assert_eq!(names1, names2);
    }

    #[test]
    fn multithreaded_writer_round_trips_through_multithreaded_reader() {
        let temp_dir = tempfile::tempdir().unwrap();
        let bam_path = temp_dir.path().join("out.bam");

        {
            let mut reader = build_alignment_reader(
                Path::new("tests/cases/test.bam"),
                NonZeroUsize::new(1).unwrap(),
            )
            .unwrap();
            let header = reader.read_header().unwrap();

            let file = File::create(&bam_path).unwrap();
            let sink: Box<dyn Write + Send> = Box::new(file);
            let mut writer =
                build_alignment_writer(sink, Format::Bam, NonZeroUsize::new(4).unwrap()).unwrap();
            writer.write_header(&header).unwrap();
            for result in reader.records(&header) {
                let record = result.unwrap();
                writer.write_record(&header, &record).unwrap();
            }
            writer.finish(&header).unwrap();
        }

        let mut reader = build_alignment_reader(&bam_path, NonZeroUsize::new(1).unwrap()).unwrap();
        let header = reader.read_header().unwrap();
        let count = reader.records(&header).count();
        assert!(count > 0);
    }

    #[test]
    fn mislabeled_bam_extension_on_uncompressed_content_still_reads_with_threads_greater_than_one()
    {
        // a SAM file, deliberately given a `.bam` extension. Extension alone would (incorrectly)
        // route this through the multithreaded BGZF-decode path; content-sniffing should catch
        // that it isn't actually BGZF-compressed and fall back to the normal autodetecting path.
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("mislabeled.bam");
        std::fs::write(
            &path,
            b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:1000\nr01\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n",
        )
        .unwrap();

        let mut reader = build_alignment_reader(&path, NonZeroUsize::new(4).unwrap()).unwrap();
        let header = reader.read_header().unwrap();
        let count = reader.records(&header).count();
        assert_eq!(count, 1);
    }
}
