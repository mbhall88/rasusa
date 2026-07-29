use crate::alignment::{is_query_grouped, program_entry};
use crate::fastx::{Fastx, FastxError, OnePassStats};
use crate::format::OutputEncoding;
use crate::subsampler::seeded_rng;
use crate::threading::build_alignment_reader;
use anyhow::Result;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::Record;
use rand::prelude::*;
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

pub trait RecordSource {
    fn read_lengths(&self) -> Result<Vec<u32>, FastxError>;
    /// Returns the number of reads (or read templates, for paired alignment data) in the source,
    /// without materializing their individual lengths.
    fn count(&self) -> Result<usize, FastxError>;
    fn filter_reads_into(
        &self,
        reads_to_keep: &[bool],
        nb_reads_keep: usize,
        write_to: &mut dyn Write,
        encoding: OutputEncoding,
    ) -> Result<usize, FastxError>;
}

pub struct AlignmentSource {
    path: PathBuf,
    threads: NonZeroUsize,
}

impl AlignmentSource {
    pub fn new(path: &Path, threads: NonZeroUsize) -> Self {
        Self {
            path: path.to_path_buf(),
            threads,
        }
    }
}

/// Extracts a record's name as an owned byte vector (empty if unset), for storage across
/// iterations. Shared by every place in this module that groups alignment records by name
/// ([`TemplateIndexer`], [`check_name_grouped`], [`TemplateGrouper`]).
fn record_name_bytes(record: &dyn Record) -> Vec<u8> {
    record
        .name()
        .map(|n| (n.as_ref() as &[u8]).to_vec())
        .unwrap_or_default()
}

/// How many records the name-grouped guard scans, when the header doesn't already declare
/// grouping/sorting, before giving up and accepting the input. This is a documented limitation
/// (see the `one_pass` field doc on [`crate::reads::Reads`]), not a guarantee: input that breaks
/// grouping only after this many records is not caught.
const NAME_GROUP_SCAN_LIMIT: usize = 50;

/// Guards one-pass SAM/BAM/CRAM streaming against silently splitting a template's records apart.
///
/// One-pass groups a template by comparing each record only to the one immediately before it -
/// correct for name-grouped input (which unaligned data off a sequencer, and collated output,
/// both are), but silently wrong on anything else. This is accepted outright when the header
/// declares grouping or name-sorting ([`is_query_grouped`]); otherwise the first
/// [`NAME_GROUP_SCAN_LIMIT`] records are scanned for positive evidence of grouping: a run of
/// consecutive records sharing a name, with no name reappearing after a different one. Input
/// with no segmented records needs no grouping and is accepted without inspecting names at all.
///
/// Only names and flags are read here - never a whole record - so peak memory doesn't grow with
/// read length or input size for SAM/BAM, whose readers decode a record's fields lazily on
/// demand. CRAM has no such lazy representation in the reader this uses (its block-based
/// compression means decoding any field forces the whole record, and often the rest of its
/// containing slice, to be decoded) - the scan is still correct for CRAM, just not
/// constant-memory in the same way.
///
/// The reader used for the scan is discarded afterwards; the caller re-opens the file for the
/// actual one-pass streaming read.
fn check_name_grouped(path: &Path, threads: NonZeroUsize) -> Result<(), FastxError> {
    let mut reader = build_alignment_reader(path, threads)
        .map_err(|source| FastxError::AlignmentReadError { source })?;
    let header = reader
        .read_header()
        .map_err(|source| FastxError::AlignmentReadError { source })?;

    if is_query_grouped(&header) {
        return Ok(());
    }

    let mut last_name: Option<Vec<u8>> = None;
    let mut closed_names: HashSet<Vec<u8>> = HashSet::new();

    for result in reader.records(&header).take(NAME_GROUP_SCAN_LIMIT) {
        let record = result.map_err(|source| FastxError::AlignmentReadError { source })?;
        let flags = record
            .flags()
            .unwrap_or(noodles::sam::alignment::record::Flags::empty());

        if !flags.is_segmented() {
            continue;
        }

        let name = record_name_bytes(&*record);

        if let Some(last) = &last_name {
            if last != &name {
                closed_names.insert(last.clone());
            }
        }
        if closed_names.contains(&name) {
            return Err(FastxError::UngroupedAlignmentInput);
        }
        last_name = Some(name);
    }

    Ok(())
}

impl AlignmentSource {
    /// Runs the name-grouped guard (see [`check_name_grouped`]) against a throwaway reader.
    ///
    /// Callers driving output through a file (rather than an in-memory sink) should run this
    /// *before* opening/truncating that file, so a rejection doesn't leave a truncated output
    /// file behind; [`Reads::run_one_pass`](crate::reads::Reads) does so.
    pub fn check_name_grouped(&self) -> Result<(), FastxError> {
        check_name_grouped(&self.path, self.threads)
    }

    /// Streams a single-pass, probabilistic fraction subsample of unaligned SAM/BAM/CRAM input:
    /// reads the input exactly once, keeping each template (a segmented read's records, grouped
    /// by comparing each record only to the one immediately before it - see
    /// [`check_name_grouped`]) independently with probability `fraction`, writing kept records
    /// immediately. Unsegmented records are each their own template.
    ///
    /// Does not itself run [`check_name_grouped`] - callers that haven't already run
    /// [`AlignmentSource::check_name_grouped`] (e.g. tests) get an unguarded stream.
    pub fn subsample_one_pass(
        &self,
        fraction: f32,
        seed: Option<u64>,
        write_to: &mut dyn Write,
        encoding: OutputEncoding,
    ) -> Result<OnePassStats, FastxError> {
        let mut reader = build_alignment_reader(&self.path, self.threads)
            .map_err(|source| FastxError::AlignmentReadError { source })?;
        let mut header = reader
            .read_header()
            .map_err(|source| FastxError::AlignmentReadError { source })?;

        let mut rng = seeded_rng(seed);
        let mut stats = OnePassStats::default();
        let mut grouper = TemplateGrouper::default();

        if let OutputEncoding::Alignment(format) = encoding {
            let (pg_id, pg_map) = program_entry(&header);
            header.programs_mut().as_mut().insert(pg_id.into(), pg_map);

            {
                let mut writer = noodles_util::alignment::io::writer::Builder::default()
                    .set_format(format)
                    .build_from_writer(&mut *write_to)
                    .map_err(|source| FastxError::AlignmentReadError { source })?;
                writer
                    .write_header(&header)
                    .map_err(|source| FastxError::AlignmentReadError { source })?;

                for result in reader.records(&header) {
                    let record =
                        result.map_err(|source| FastxError::AlignmentReadError { source })?;
                    let (is_new_template, decision) =
                        grouper.decide(&*record, &mut rng, fraction)?;
                    if is_new_template {
                        stats.reads_seen += 1;
                        if decision {
                            stats.reads_kept += 1;
                        }
                    }
                    if decision {
                        writer
                            .write_record(&header, &record)
                            .map_err(|source| FastxError::AlignmentReadError { source })?;
                    }
                }
                writer
                    .finish(&header)
                    .map_err(|source| FastxError::AlignmentReadError { source })?;
            }
            write_to.flush().map_err(|source| FastxError::WriteError {
                source: anyhow::Error::from(source),
            })?;
        } else {
            let is_fasta = matches!(encoding, OutputEncoding::Fastx { fasta: true });

            for result in reader.records(&header) {
                let record = result.map_err(|source| FastxError::AlignmentReadError { source })?;
                let (is_new_template, decision) = grouper.decide(&*record, &mut rng, fraction)?;
                if is_new_template {
                    stats.reads_seen += 1;
                    if decision {
                        stats.reads_kept += 1;
                    }
                }
                if decision {
                    let name = record.name().map(|n| n.as_ref()).unwrap_or(&b"*"[..]);
                    let seq: Vec<u8> = record.sequence().iter().collect();
                    let qual: Vec<u8> = record
                        .quality_scores()
                        .iter()
                        .map(|q| q.map(|score| score + 33))
                        .collect::<Result<Vec<u8>, _>>()
                        .map_err(|source| FastxError::AlignmentReadError { source })?;
                    let qual = if qual.is_empty() {
                        None
                    } else {
                        Some(&qual[..])
                    };
                    crate::record::write_fastx_record(write_to, name, &seq, qual, is_fasta, b"\n")?;
                }
            }
        }

        Ok(stats)
    }
}

/// Groups consecutive alignment records into templates by comparing each record only to the one
/// immediately before it - the same constant-memory grouping [`check_name_grouped`] guards - and
/// carries one keep/drop decision per template forward across the records that share it.
#[derive(Default)]
struct TemplateGrouper {
    current_name: Option<Vec<u8>>,
    current_decision: bool,
}

impl TemplateGrouper {
    /// Decides whether `record` belongs to a new template (relative to the previous call's
    /// record) and, if so, draws a fresh keep/drop decision; otherwise reuses the previous
    /// decision so every record of a template is kept or dropped together. Also rejects mapped
    /// reads, matching existing (two-pass) behaviour.
    ///
    /// Returns `(is_new_template, decision)`.
    fn decide(
        &mut self,
        record: &dyn Record,
        rng: &mut impl Rng,
        fraction: f32,
    ) -> Result<(bool, bool), FastxError> {
        let flags = record
            .flags()
            .map_err(|source| FastxError::AlignmentReadError { source })?;
        if !flags.is_unmapped() {
            return Err(FastxError::MappedReadDetected);
        }

        let name = flags.is_segmented().then(|| record_name_bytes(record));

        let is_new_template = match (&name, self.current_name.as_ref()) {
            (Some(n), Some(cur)) => n != cur,
            _ => true,
        };

        if is_new_template {
            self.current_decision = rng.random_bool(fraction as f64);
        }
        self.current_name = name;

        Ok((is_new_template, self.current_decision))
    }
}

/// Assigns a stable per-template index to alignment records, coalescing segmented (paired)
/// records that share a QNAME onto the same index the way SAM template grouping expects.
/// Indices are handed out sequentially (`0, 1, 2, ...`) in first-seen order, so `next_idx` at any
/// point also equals the number of distinct templates seen so far.
///
/// Every alignment-format scan in this module (`read_lengths`, `count`, `filter_reads_into`)
/// needs exactly this grouping, so it's centralised here rather than reimplemented per scan.
#[derive(Default)]
struct TemplateIndexer {
    next_idx: usize,
    qname_to_idx: HashMap<Vec<u8>, usize>,
}

impl TemplateIndexer {
    fn index_for(&mut self, record: &dyn Record, flags: Flags) -> usize {
        if !flags.is_segmented() {
            let idx = self.next_idx;
            self.next_idx += 1;
            return idx;
        }

        let name = record_name_bytes(record);

        if let Some(&idx) = self.qname_to_idx.get(&name) {
            return idx;
        }
        let idx = self.next_idx;
        self.next_idx += 1;
        self.qname_to_idx.insert(name, idx);
        idx
    }
}

impl RecordSource for AlignmentSource {
    fn read_lengths(&self) -> Result<Vec<u32>, FastxError> {
        let mut reader = build_alignment_reader(&self.path, self.threads)
            .map_err(|source| FastxError::AlignmentReadError { source })?;

        let header = reader
            .read_header()
            .map_err(|source| FastxError::AlignmentReadError { source })?;

        let mut read_lengths: Vec<u32> = vec![];
        let mut indexer = TemplateIndexer::default();

        for result in reader.records(&header) {
            let record = result.map_err(|source| FastxError::AlignmentReadError { source })?;
            let flags = record
                .flags()
                .unwrap_or(noodles::sam::alignment::record::Flags::empty());

            if !flags.is_unmapped() {
                return Err(FastxError::MappedReadDetected);
            }

            let rlen = record.sequence().len() as u32;
            let idx = indexer.index_for(&*record, flags);
            if idx == read_lengths.len() {
                read_lengths.push(rlen);
            } else {
                read_lengths[idx] += rlen;
            }
        }

        Ok(read_lengths)
    }

    fn count(&self) -> Result<usize, FastxError> {
        let mut reader = build_alignment_reader(&self.path, self.threads)
            .map_err(|source| FastxError::AlignmentReadError { source })?;

        let header = reader
            .read_header()
            .map_err(|source| FastxError::AlignmentReadError { source })?;

        let mut indexer = TemplateIndexer::default();

        for result in reader.records(&header) {
            let record = result.map_err(|source| FastxError::AlignmentReadError { source })?;
            let flags = record
                .flags()
                .unwrap_or(noodles::sam::alignment::record::Flags::empty());

            if !flags.is_unmapped() {
                return Err(FastxError::MappedReadDetected);
            }

            indexer.index_for(&*record, flags);
        }

        // indices are assigned sequentially in first-seen order, so the number of templates seen
        // equals the count of indices handed out.
        Ok(indexer.next_idx)
    }

    fn filter_reads_into(
        &self,
        reads_to_keep: &[bool],
        nb_reads_keep: usize,
        write_to: &mut dyn Write,
        encoding: OutputEncoding,
    ) -> Result<usize, FastxError> {
        let mut reader = build_alignment_reader(&self.path, self.threads)
            .map_err(|source| FastxError::AlignmentReadError { source })?;

        let mut header = reader
            .read_header()
            .map_err(|source| FastxError::AlignmentReadError { source })?;

        let mut total_len = 0;
        let mut written_templates: HashSet<usize> = HashSet::new();
        let mut indexer = TemplateIndexer::default();

        // If the output format is an alignment format (SAM/BAM/CRAM), we write alignment records
        if let OutputEncoding::Alignment(format) = encoding {
            {
                // add rasusa program command line to header
                let (pg_id, pg_map) = program_entry(&header);
                header.programs_mut().as_mut().insert(pg_id.into(), pg_map);

                // `write_to` is a borrowed `&mut dyn Write`, not an owned `Send` sink, so it can't
                // be wrapped in the multithreaded BGZF encoder from `crate::threading` (that
                // needs to own the writer for its background compression threads). BAM output
                // written through `reads` is therefore always single-threaded; `--threads` only
                // speeds up BAM *reading* here.
                let mut writer = noodles_util::alignment::io::writer::Builder::default()
                    .set_format(format)
                    .build_from_writer(&mut *write_to)
                    .map_err(|source| FastxError::AlignmentReadError { source })?;

                writer
                    .write_header(&header)
                    .map_err(|source| FastxError::AlignmentReadError { source })?;

                for result in reader.records(&header) {
                    let record =
                        result.map_err(|source| FastxError::AlignmentReadError { source })?;
                    let flags = record
                        .flags()
                        .unwrap_or(noodles::sam::alignment::record::Flags::empty());
                    let current_idx = indexer.index_for(&*record, flags);

                    if current_idx < reads_to_keep.len() && reads_to_keep[current_idx] {
                        total_len += record.sequence().len();
                        writer
                            .write_record(&header, &record)
                            .map_err(|source| FastxError::AlignmentReadError { source })?;
                        written_templates.insert(current_idx);
                    }
                }
                writer
                    .finish(&header)
                    .map_err(|source| FastxError::AlignmentReadError { source })?;
            }
            write_to.flush().map_err(|source| FastxError::WriteError {
                source: anyhow::Error::from(source),
            })?;
        } else {
            let is_fasta = matches!(encoding, OutputEncoding::Fastx { fasta: true });
            // Otherwise, we output as FASTQ (or FASTA)
            for result in reader.records(&header) {
                let record = result.map_err(|source| FastxError::AlignmentReadError { source })?;
                let flags = record
                    .flags()
                    .unwrap_or(noodles::sam::alignment::record::Flags::empty());
                let current_idx = indexer.index_for(&*record, flags);

                if current_idx < reads_to_keep.len() && reads_to_keep[current_idx] {
                    total_len += record.sequence().len();

                    let name = record.name().map(|n| n.as_ref()).unwrap_or(&b"*"[..]);
                    let seq: Vec<u8> = record.sequence().iter().collect();
                    let qual: Vec<u8> = record
                        .quality_scores()
                        .iter()
                        .map(|q| q.map(|score| score + 33))
                        .collect::<Result<Vec<u8>, _>>()
                        .map_err(|source| FastxError::AlignmentReadError { source })?;
                    let qual = if qual.is_empty() {
                        None
                    } else {
                        Some(&qual[..])
                    };

                    crate::record::write_fastx_record(write_to, name, &seq, qual, is_fasta, b"\n")?;

                    written_templates.insert(current_idx);
                }
            }
        }

        if written_templates.len() == nb_reads_keep {
            Ok(total_len)
        } else {
            Err(FastxError::IndicesNotFound)
        }
    }
}

pub fn determine_record_source(path: &Path, threads: NonZeroUsize) -> Box<dyn RecordSource> {
    match path.extension().and_then(|ext| ext.to_str()) {
        Some("sam") | Some("bam") | Some("cram") => Box::new(AlignmentSource::new(path, threads)),
        _ => Box::new(Fastx::from_path(path)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_util::alignment::io::Format;
    use std::fs::File;
    use std::path::Path;

    fn create_test_sam(path: &Path) {
        let mut file = File::create(path).unwrap();
        let content = b"@HD\tVN:1.6\tSO:coordinate\n\
                        @SQ\tSN:ref\tLN:1000\n\
                        r01\t4\tref\t10\t255\t10M\t*\t0\t0\tACTGACTGAW\t*\n\
                        r02\t4\tref\t20\t255\t10M\t*\t0\t0\tGCTGACTGAC\t*\n";
        file.write_all(content).unwrap();
    }

    fn create_test_paired_sam(path: &Path) {
        let mut file = File::create(path).unwrap();
        let content = b"@HD\tVN:1.6\tSO:coordinate\n\
                        @SQ\tSN:ref\tLN:1000\n\
                        r01\t77\tref\t10\t255\t10M\t=\t20\t20\tACTGACTGAC\t*\n\
                        r01\t141\tref\t20\t255\t10M\t=\t10\t-20\tGCTGACTGAC\t*\n\
                        r02\t77\tref\t30\t255\t10M\t=\t40\t20\tACTGACTGAC\t*\n\
                        r02\t141\tref\t40\t255\t10M\t=\t30\t-20\tGCTGACTGAC\t*\n";
        file.write_all(content).unwrap();
    }

    fn create_test_sam_with_quality(path: &Path) {
        let mut file = File::create(path).unwrap();
        let content = b"@HD\tVN:1.6\tSO:coordinate\n\
                        @SQ\tSN:ref\tLN:1000\n\
                        r01\t4\tref\t10\t255\t10M\t*\t0\t0\tACTGACTGAC\tIIIIIIIIII\n\
                        r02\t4\tref\t20\t255\t10M\t*\t0\t0\tGCTGACTGAC\tJJJJJJJJJJ\n";
        file.write_all(content).unwrap();
    }

    #[test]
    fn test_determine_record_source_alignment() {
        let path = Path::new("test.bam");
        let _source = determine_record_source(path, NonZeroUsize::new(1).unwrap());
    }

    #[test]
    fn test_alignment_source_read_lengths() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("test.sam");
        create_test_sam(&path);

        let source = AlignmentSource::new(&path, NonZeroUsize::new(1).unwrap());
        let actual = source.read_lengths().unwrap();
        assert_eq!(actual.len(), 2);
        assert_eq!(actual[0], 10);
        assert_eq!(actual[1], 10);
    }

    #[test]
    fn test_alignment_source_read_lengths_paired() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("test_paired.sam");
        create_test_paired_sam(&path);

        let source = AlignmentSource::new(&path, NonZeroUsize::new(1).unwrap());
        let actual = source.read_lengths().unwrap();
        // 2 templates, each length 10 + 10 = 20.
        assert_eq!(actual.len(), 2);
        assert_eq!(actual[0], 20);
        assert_eq!(actual[1], 20);
    }

    #[test]
    fn test_alignment_source_count() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("test.sam");
        create_test_sam(&path);

        let source = AlignmentSource::new(&path, NonZeroUsize::new(1).unwrap());
        assert_eq!(source.count().unwrap(), 2);
    }

    #[test]
    fn test_alignment_source_count_paired_dedups_by_template() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("test_paired.sam");
        create_test_paired_sam(&path);

        let source = AlignmentSource::new(&path, NonZeroUsize::new(1).unwrap());
        // 4 records but 2 templates (read pairs)
        assert_eq!(source.count().unwrap(), 2);
    }

    #[test]
    fn test_alignment_source_filter_reads_into_sam_to_bam() {
        let temp_dir = tempfile::tempdir().unwrap();
        let input_path = temp_dir.path().join("test.sam");
        create_test_sam(&input_path);

        let source = AlignmentSource::new(&input_path, NonZeroUsize::new(1).unwrap());
        let reads_to_keep = vec![true, false];
        let nb_reads_keep = 1;
        let mut buffer = Vec::new();

        let result = source.filter_reads_into(
            &reads_to_keep,
            nb_reads_keep,
            &mut buffer,
            OutputEncoding::Alignment(Format::Bam),
        );
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), 10);

        // Verify the output is a valid BAM and has 1 record
        let mut reader = noodles_util::alignment::io::reader::Builder::default()
            .build_from_reader(&buffer[..])
            .unwrap();
        let header = reader.read_header().unwrap();
        let mut count = 0;
        for result in reader.records(&header) {
            let record = result.unwrap();
            assert_eq!(record.name().map(|n| n.as_ref()), Some(&b"r01"[..]));
            count += 1;
        }
        assert_eq!(count, 1);
    }

    #[test]
    fn test_alignment_source_filter_reads_into_sam_to_fasta() {
        let temp_dir = tempfile::tempdir().unwrap();
        let input_path = temp_dir.path().join("test.sam");
        create_test_sam(&input_path);

        let source = AlignmentSource::new(&input_path, NonZeroUsize::new(1).unwrap());
        let reads_to_keep = vec![false, true];
        let nb_reads_keep = 1;
        let mut buffer = Vec::new();

        let result = source.filter_reads_into(
            &reads_to_keep,
            nb_reads_keep,
            &mut buffer,
            OutputEncoding::Fastx { fasta: true },
        );
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), 10);

        let output = String::from_utf8(buffer).unwrap();
        assert_eq!(output.lines().next().unwrap(), ">r02");
        assert_eq!(output.lines().count(), 2);
        assert!(output.contains("GCTGACTGAC"));
    }

    #[test]
    fn test_alignment_source_filter_reads_into_sam_to_fastq_with_quality() {
        let temp_dir = tempfile::tempdir().unwrap();
        let input_path = temp_dir.path().join("test_qual.sam");
        create_test_sam_with_quality(&input_path);

        let source = AlignmentSource::new(&input_path, NonZeroUsize::new(1).unwrap());
        let reads_to_keep = vec![false, true];
        let nb_reads_keep = 1;
        let mut buffer = Vec::new();

        let result = source.filter_reads_into(
            &reads_to_keep,
            nb_reads_keep,
            &mut buffer,
            OutputEncoding::Fastx { fasta: false },
        );
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), 10);

        let output = String::from_utf8(buffer).unwrap();
        let lines: Vec<&str> = output.lines().collect();
        assert_eq!(lines[0], "@r02");
        assert_eq!(lines[1], "GCTGACTGAC");
        assert_eq!(lines[2], "+");
        assert_eq!(lines[3], "JJJJJJJJJJ");
    }

    #[test]
    fn test_alignment_source_filter_reads_into_paired_sam_to_sam() {
        let temp_dir = tempfile::tempdir().unwrap();
        let input_path = temp_dir.path().join("test_paired.sam");
        create_test_paired_sam(&input_path);

        let source = AlignmentSource::new(&input_path, NonZeroUsize::new(1).unwrap());
        // 2 templates. Keep first.
        let reads_to_keep = vec![true, false];
        let nb_reads_keep = 1;
        let mut buffer = Vec::new();

        let result = source.filter_reads_into(
            &reads_to_keep,
            nb_reads_keep,
            &mut buffer,
            OutputEncoding::Alignment(Format::Sam),
        );
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), 20);

        // Verify the output has 2 records (1 pair)
        let mut reader = noodles_util::alignment::io::reader::Builder::default()
            .build_from_reader(&buffer[..])
            .unwrap();
        let header = reader.read_header().unwrap();
        let mut count = 0;
        for result in reader.records(&header) {
            let record = result.unwrap();
            assert_eq!(record.name().map(|n| n.as_ref()), Some(&b"r01"[..]));
            count += 1;
        }
        assert_eq!(count, 2);
    }

    #[test]
    fn test_alignment_source_filter_reads_into_bam_to_sam() {
        let temp_dir = tempfile::tempdir().unwrap();
        let sam_path = temp_dir.path().join("input.sam");
        create_test_sam_with_quality(&sam_path);

        // Convert SAM to BAM first
        let bam_path = temp_dir.path().join("input.bam");
        {
            let source = AlignmentSource::new(&sam_path, NonZeroUsize::new(1).unwrap());
            let mut bam_file = File::create(&bam_path).unwrap();
            source
                .filter_reads_into(
                    &[true, true],
                    2,
                    &mut bam_file,
                    OutputEncoding::Alignment(Format::Bam),
                )
                .unwrap();
        }

        let source = AlignmentSource::new(&bam_path, NonZeroUsize::new(1).unwrap());
        let reads_to_keep = vec![false, true];
        let nb_reads_keep = 1;
        let mut buffer = Vec::new();

        let result = source.filter_reads_into(
            &reads_to_keep,
            nb_reads_keep,
            &mut buffer,
            OutputEncoding::Alignment(Format::Sam),
        );
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), 10);

        let output = String::from_utf8(buffer).unwrap();
        assert!(output.contains("r02"));
        assert!(output.contains("GCTGACTGAC"));
        assert!(output.contains("JJJJJJJJJJ"));
    }

    fn create_test_sam_with_header(path: &Path, header_extra: &str) {
        let mut file = File::create(path).unwrap();
        let content = format!(
            "@HD\tVN:1.6\t{header_extra}\n\
             @SQ\tSN:ref\tLN:1000\n\
             r01\t77\tref\t10\t255\t10M\t=\t20\t20\tACTGACTGAC\t*\n\
             r02\t77\tref\t30\t255\t10M\t=\t40\t20\tACTGACTGAC\t*\n\
             r01\t141\tref\t20\t255\t10M\t=\t10\t-20\tGCTGACTGAC\t*\n\
             r02\t141\tref\t40\t255\t10M\t=\t30\t-20\tGCTGACTGAC\t*\n"
        );
        file.write_all(content.as_bytes()).unwrap();
    }

    #[test]
    fn check_name_grouped_accepts_adjacent_pairs() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("grouped.sam");
        create_test_paired_sam(&path);

        assert!(check_name_grouped(&path, NonZeroUsize::new(1).unwrap()).is_ok());
    }

    #[test]
    fn check_name_grouped_accepts_unsegmented_records() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("unsegmented.sam");
        create_test_sam(&path);

        assert!(check_name_grouped(&path, NonZeroUsize::new(1).unwrap()).is_ok());
    }

    #[test]
    fn check_name_grouped_rejects_non_adjacent_pairs() {
        // r01's two records are split apart by r02's - not name-grouped.
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("non_adjacent.sam");
        create_test_sam_with_header(&path, "SO:coordinate");

        let err = check_name_grouped(&path, NonZeroUsize::new(1).unwrap())
            .expect_err("non-adjacent mates should be rejected");
        assert!(matches!(err, FastxError::UngroupedAlignmentInput));
    }

    #[test]
    fn check_name_grouped_trusts_go_query_header_without_scanning() {
        // Same non-adjacent layout that `check_name_grouped_rejects_non_adjacent_pairs` rejects,
        // but a `GO:query` header short-circuits the scan entirely.
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("trusted.sam");
        create_test_sam_with_header(&path, "GO:query");

        assert!(check_name_grouped(&path, NonZeroUsize::new(1).unwrap()).is_ok());
    }

    #[test]
    fn check_name_grouped_trusts_so_queryname_header_without_scanning() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("trusted.sam");
        create_test_sam_with_header(&path, "SO:queryname");

        assert!(check_name_grouped(&path, NonZeroUsize::new(1).unwrap()).is_ok());
    }

    #[test]
    fn test_alignment_source_subsample_one_pass_fraction_one_keeps_every_record() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("test.sam");
        create_test_sam(&path);

        let source = AlignmentSource::new(&path, NonZeroUsize::new(1).unwrap());
        let mut out = Vec::new();
        let stats = source
            .subsample_one_pass(
                1.0,
                Some(1),
                &mut out,
                OutputEncoding::Fastx { fasta: false },
            )
            .unwrap();

        assert_eq!(stats.reads_seen, 2);
        assert_eq!(stats.reads_kept, 2);
        let output = String::from_utf8(out).unwrap();
        assert!(output.contains("r01"));
        assert!(output.contains("r02"));
    }

    #[test]
    fn test_alignment_source_subsample_one_pass_fraction_zero_keeps_nothing() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("test.sam");
        create_test_sam(&path);

        let source = AlignmentSource::new(&path, NonZeroUsize::new(1).unwrap());
        let mut out = Vec::new();
        let stats = source
            .subsample_one_pass(
                0.0,
                Some(1),
                &mut out,
                OutputEncoding::Fastx { fasta: false },
            )
            .unwrap();

        assert_eq!(stats.reads_seen, 2);
        assert_eq!(stats.reads_kept, 0);
        assert!(out.is_empty());
    }

    #[test]
    fn test_alignment_source_subsample_one_pass_keeps_mates_together() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("test_paired.sam");
        create_test_paired_sam(&path);

        let source = AlignmentSource::new(&path, NonZeroUsize::new(1).unwrap());
        let mut out = Vec::new();
        let stats = source
            .subsample_one_pass(
                1.0,
                Some(1),
                &mut out,
                OutputEncoding::Fastx { fasta: false },
            )
            .unwrap();

        // 2 templates (r01, r02), each contributing 2 records.
        assert_eq!(stats.reads_seen, 2);
        assert_eq!(stats.reads_kept, 2);
        let output = String::from_utf8(out).unwrap();
        assert_eq!(output.matches("r01").count(), 2);
        assert_eq!(output.matches("r02").count(), 2);
    }

    #[test]
    fn test_alignment_source_subsample_one_pass_rejects_mapped_reads() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("mapped.sam");
        let mut file = File::create(&path).unwrap();
        file.write_all(
            b"@HD\tVN:1.6\tSO:coordinate\n\
              @SQ\tSN:ref\tLN:1000\n\
              r01\t0\tref\t10\t255\t10M\t*\t0\t0\tACTGACTGAC\t*\n",
        )
        .unwrap();

        let source = AlignmentSource::new(&path, NonZeroUsize::new(1).unwrap());
        let mut out = Vec::new();
        let result = source.subsample_one_pass(
            1.0,
            Some(1),
            &mut out,
            OutputEncoding::Fastx { fasta: false },
        );

        assert!(matches!(result, Err(FastxError::MappedReadDetected)));
    }

    #[test]
    fn test_alignment_source_check_name_grouped_rejects_ungrouped_input() {
        // subsample_one_pass doesn't run the guard itself (see its doc comment) - callers run
        // AlignmentSource::check_name_grouped first, as Reads::run_one_pass does, so a rejection
        // happens before the output file is created/truncated.
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("non_adjacent.sam");
        create_test_sam_with_header(&path, "SO:coordinate");

        let source = AlignmentSource::new(&path, NonZeroUsize::new(1).unwrap());
        let result = source.check_name_grouped();

        assert!(matches!(result, Err(FastxError::UngroupedAlignmentInput)));
    }

    #[test]
    fn test_alignment_source_subsample_one_pass_same_seed_gives_same_output() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("test_paired.sam");
        create_test_paired_sam(&path);

        let source = AlignmentSource::new(&path, NonZeroUsize::new(1).unwrap());

        let mut out1 = Vec::new();
        source
            .subsample_one_pass(
                0.5,
                Some(7),
                &mut out1,
                OutputEncoding::Fastx { fasta: false },
            )
            .unwrap();

        let mut out2 = Vec::new();
        source
            .subsample_one_pass(
                0.5,
                Some(7),
                &mut out2,
                OutputEncoding::Fastx { fasta: false },
            )
            .unwrap();

        assert_eq!(out1, out2);
    }

    #[test]
    fn test_alignment_source_subsample_one_pass_alignment_output_round_trips() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("test_paired.sam");
        create_test_paired_sam(&path);

        let source = AlignmentSource::new(&path, NonZeroUsize::new(1).unwrap());
        let mut out = Vec::new();
        let stats = source
            .subsample_one_pass(
                1.0,
                Some(1),
                &mut out,
                OutputEncoding::Alignment(Format::Sam),
            )
            .unwrap();
        assert_eq!(stats.reads_kept, 2);

        let mut reader = noodles_util::alignment::io::reader::Builder::default()
            .build_from_reader(&out[..])
            .unwrap();
        let header = reader.read_header().unwrap();
        let count = reader.records(&header).count();
        assert_eq!(count, 4);
    }
}
