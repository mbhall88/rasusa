use crate::fastx::{Fastx, FastxError};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::collections::{HashMap, HashSet};
use noodles_util::alignment::io::Format;
use anyhow::Result;

pub trait RecordSource {
    fn read_lengths(&self) -> Result<Vec<u32>, FastxError>;
    fn filter_reads_into(
        &self,
        reads_to_keep: &[bool],
        nb_reads_keep: usize,
        write_to: &mut dyn Write,
        output_format: Option<Format>,
    ) -> Result<usize, FastxError>;
}

pub struct AlignmentSource {
    path: PathBuf,
}

impl AlignmentSource {
    pub fn new(path: &Path) -> Self {
        Self {
            path: path.to_path_buf(),
        }
    }
}

impl RecordSource for AlignmentSource {
    fn read_lengths(&self) -> Result<Vec<u32>, FastxError> {
        let mut reader = noodles_util::alignment::io::reader::Builder::default()
            .build_from_path(&self.path)
            .map_err(|source| FastxError::AlignmentReadError { source })?;

        let header = reader.read_header()
            .map_err(|source| FastxError::AlignmentReadError { source })?;

        let mut read_lengths: Vec<u32> = vec![];
        let mut qname_to_idx: HashMap<Vec<u8>, usize> = HashMap::new();

        for result in reader.records(&header) {
            let record = result.map_err(|source| FastxError::AlignmentReadError { source })?;
            let flags = record.flags().unwrap_or(noodles::sam::alignment::record::Flags::empty());
            
            if flags.is_segmented() {
                let name = record.name().map(|n| (n.as_ref() as &[u8]).to_vec()).unwrap_or_default();
                if let Some(&idx) = qname_to_idx.get(&name) {
                    read_lengths[idx] += record.sequence().len() as u32;
                } else {
                    qname_to_idx.insert(name, read_lengths.len());
                    read_lengths.push(record.sequence().len() as u32);
                }
            } else {
                read_lengths.push(record.sequence().len() as u32);
            }
        }

        Ok(read_lengths)
    }

    fn filter_reads_into(
        &self,
        reads_to_keep: &[bool],
        nb_reads_keep: usize,
        write_to: &mut dyn Write,
        output_format: Option<Format>,
    ) -> Result<usize, FastxError> {
        let mut reader = noodles_util::alignment::io::reader::Builder::default()
            .build_from_path(&self.path)
            .map_err(|source| FastxError::AlignmentReadError { source })?;

        let header = reader.read_header()
            .map_err(|source| FastxError::AlignmentReadError { source })?;

        let mut read_idx: usize = 0;
        let mut total_len = 0;
        let mut written_templates: HashSet<usize> = HashSet::new();
        let mut qname_to_idx: HashMap<Vec<u8>, usize> = HashMap::new();

        // If the output format is an alignment format (SAM/BAM/CRAM), we write alignment records
        if let Some(format) = output_format {
            let mut writer = noodles_util::alignment::io::writer::Builder::default()
                .set_format(format)
                .build_from_writer(write_to)
                .map_err(|source| FastxError::AlignmentReadError { source })?;

            writer.write_header(&header)
                .map_err(|source| FastxError::AlignmentReadError { source })?;

            for result in reader.records(&header) {
                let record = result.map_err(|source| FastxError::AlignmentReadError { source })?;
                let flags = record.flags().unwrap_or(noodles::sam::alignment::record::Flags::empty());

                let current_idx = if flags.is_segmented() {
                    let name = record.name().map(|n| (n.as_ref() as &[u8]).to_vec()).unwrap_or_default();
                    if let Some(&idx) = qname_to_idx.get(&name) {
                        idx
                    } else {
                        let idx = read_idx;
                        qname_to_idx.insert(name, idx);
                        read_idx += 1;
                        idx
                    }
                } else {
                    let idx = read_idx;
                    read_idx += 1;
                    idx
                };

                if current_idx < reads_to_keep.len() && reads_to_keep[current_idx] {
                    total_len += record.sequence().len();
                    writer.write_record(&header, &record)
                        .map_err(|source| FastxError::AlignmentReadError { source })?;
                    written_templates.insert(current_idx);
                }
            }
        } else {
            // Otherwise, we output as FASTQ (or FASTA)
            for result in reader.records(&header) {
                let record = result.map_err(|source| FastxError::AlignmentReadError { source })?;
                let flags = record.flags().unwrap_or(noodles::sam::alignment::record::Flags::empty());

                let current_idx = if flags.is_segmented() {
                    let name = record.name().map(|n| (n.as_ref() as &[u8]).to_vec()).unwrap_or_default();
                    if let Some(&idx) = qname_to_idx.get(&name) {
                        idx
                    } else {
                        let idx = read_idx;
                        qname_to_idx.insert(name, idx);
                        read_idx += 1;
                        idx
                    }
                } else {
                    let idx = read_idx;
                    read_idx += 1;
                    idx
                };

                if current_idx < reads_to_keep.len() && reads_to_keep[current_idx] {
                    total_len += record.sequence().len();
                    
                    let name = record.name().map(|n| n.as_ref()).unwrap_or(&b"*"[..]);
                    let seq: Vec<u8> = record.sequence().iter().collect();
                    let qual: Vec<u8> = record.quality_scores().iter().map(|q| q.map(|score| score + 33)).collect::<Result<Vec<u8>, _>>()
                        .map_err(|source| FastxError::AlignmentReadError { source })?;

                    if qual.is_empty() {
                        // FASTA
                        write_to.write_all(b">")
                            .map_err(|source| FastxError::WriteError { source: anyhow::Error::from(source) })?;
                        write_to.write_all(name)
                            .map_err(|source| FastxError::WriteError { source: anyhow::Error::from(source) })?;
                        write_to.write_all(b"\n")
                            .map_err(|source| FastxError::WriteError { source: anyhow::Error::from(source) })?;
                        write_to.write_all(&seq)
                            .map_err(|source| FastxError::WriteError { source: anyhow::Error::from(source) })?;
                        write_to.write_all(b"\n")
                            .map_err(|source| FastxError::WriteError { source: anyhow::Error::from(source) })?;
                    } else {
                        // FASTQ
                        write_to.write_all(b"@")
                            .map_err(|source| FastxError::WriteError { source: anyhow::Error::from(source) })?;
                        write_to.write_all(name)
                            .map_err(|source| FastxError::WriteError { source: anyhow::Error::from(source) })?;
                        write_to.write_all(b"\n")
                            .map_err(|source| FastxError::WriteError { source: anyhow::Error::from(source) })?;
                        write_to.write_all(&seq)
                            .map_err(|source| FastxError::WriteError { source: anyhow::Error::from(source) })?;
                        write_to.write_all(b"\n+\n")
                            .map_err(|source| FastxError::WriteError { source: anyhow::Error::from(source) })?;
                        write_to.write_all(&qual)
                            .map_err(|source| FastxError::WriteError { source: anyhow::Error::from(source) })?;
                        write_to.write_all(b"\n")
                            .map_err(|source| FastxError::WriteError { source: anyhow::Error::from(source) })?;
                    }

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

pub fn determine_record_source(path: &Path) -> Box<dyn RecordSource> {
    match path.extension().and_then(|ext| ext.to_str()) {
        Some("sam") | Some("bam") | Some("cram") => Box::new(AlignmentSource::new(path)),
        _ => Box::new(Fastx::from_path(path)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn test_determine_record_source_alignment() {
        let path = Path::new("test.bam");
        let _source = determine_record_source(path);
    }

    #[test]
    fn test_alignment_source_read_lengths() {
        let path = Path::new("tests/cases/test.bam");
        let source = AlignmentSource::new(path);
        let actual = source.read_lengths().unwrap();
        assert!(!actual.is_empty());
    }
}
