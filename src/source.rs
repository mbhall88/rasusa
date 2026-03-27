use crate::fastx::{Fastx, FastxError};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::collections::{HashMap, HashSet};
use noodles_util::alignment::io::Format;
use anyhow::Result;
use noodles::sam::header::record::value::map::{Program, program::tag};
use noodles::sam::header::record::value::Map;
use crate::alignment::make_program_id_unique;

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

    fn program_entry(&self, header: &noodles::sam::Header) -> (String, Map<Program>) {
        let (program_id, previous_pgid) = make_program_id_unique(header, "rasusa");

        let mut record = Map::<Program>::builder();
        record = record.insert(tag::NAME, "rasusa");
        record = record.insert(tag::VERSION, env!("CARGO_PKG_VERSION"));

        let cl = std::env::args().collect::<Vec<String>>().join(" ");
        record = record.insert(tag::COMMAND_LINE, cl);

        if let Some(pp) = previous_pgid {
            record = record.insert(tag::PREVIOUS_PROGRAM_ID, pp);
        };

        let program = record.build().expect("Failed to build program record");
        (program_id.into_owned(), program)
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
            
            if !flags.is_unmapped() {
                return Err(FastxError::MappedReadDetected);
            }

            let rlen = record.sequence().len() as u32;
            
            if flags.is_segmented() {
                let name = record.name().map(|n| (n.as_ref() as &[u8]).to_vec()).unwrap_or_default();
                if let Some(&idx) = qname_to_idx.get(&name) {
                    read_lengths[idx] += rlen;
                } else {
                    qname_to_idx.insert(name, read_lengths.len());
                    read_lengths.push(rlen);
                }
            } else {
                read_lengths.push(rlen);
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

        let mut header = reader.read_header()
            .map_err(|source| FastxError::AlignmentReadError { source })?;

        let mut read_idx: usize = 0;
        let mut total_len = 0;
        let mut written_templates: HashSet<usize> = HashSet::new();
        let mut qname_to_idx: HashMap<Vec<u8>, usize> = HashMap::new();

        // If the output format is an alignment format (SAM/BAM/CRAM), we write alignment records
        if let Some(format) = output_format {
            {
                // add rasusa program command line to header
                let (pg_id, pg_map) = self.program_entry(&header);
                header.programs_mut().as_mut().insert(pg_id.into(), pg_map);

                let mut writer = noodles_util::alignment::io::writer::Builder::default()
                    .set_format(format)
                    .build_from_writer(&mut *write_to)
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
                writer.finish(&header)
                    .map_err(|source| FastxError::AlignmentReadError { source })?;
            }
            write_to.flush()
                .map_err(|source| FastxError::WriteError { source: anyhow::Error::from(source) })?;
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
    use std::fs::File;
    use noodles_util::alignment::io::Format;

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
        let _source = determine_record_source(path);
    }

    #[test]
    fn test_alignment_source_read_lengths() {
        let temp_dir = tempfile::tempdir().unwrap();
        let path = temp_dir.path().join("test.sam");
        create_test_sam(&path);
        
        let source = AlignmentSource::new(&path);
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

        let source = AlignmentSource::new(&path);
        let actual = source.read_lengths().unwrap();
        // 2 templates, each length 10 + 10 = 20.
        assert_eq!(actual.len(), 2);
        assert_eq!(actual[0], 20);
        assert_eq!(actual[1], 20);
    }

    #[test]
    fn test_alignment_source_filter_reads_into_sam_to_bam() {
        let temp_dir = tempfile::tempdir().unwrap();
        let input_path = temp_dir.path().join("test.sam");
        create_test_sam(&input_path);
        
        let source = AlignmentSource::new(&input_path);
        let reads_to_keep = vec![true, false];
        let nb_reads_keep = 1;
        let mut buffer = Vec::new();
        
        let result = source.filter_reads_into(&reads_to_keep, nb_reads_keep, &mut buffer, Some(Format::Bam));
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
        
        let source = AlignmentSource::new(&input_path);
        let reads_to_keep = vec![false, true];
        let nb_reads_keep = 1;
        let mut buffer = Vec::new();
        
        let result = source.filter_reads_into(&reads_to_keep, nb_reads_keep, &mut buffer, None);
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
        
        let source = AlignmentSource::new(&input_path);
        let reads_to_keep = vec![false, true];
        let nb_reads_keep = 1;
        let mut buffer = Vec::new();
        
        let result = source.filter_reads_into(&reads_to_keep, nb_reads_keep, &mut buffer, None);
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
        
        let source = AlignmentSource::new(&input_path);
        // 2 templates. Keep first.
        let reads_to_keep = vec![true, false];
        let nb_reads_keep = 1;
        let mut buffer = Vec::new();
        
        let result = source.filter_reads_into(&reads_to_keep, nb_reads_keep, &mut buffer, Some(Format::Sam));
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
            let source = AlignmentSource::new(&sam_path);
            let mut bam_file = File::create(&bam_path).unwrap();
            source.filter_reads_into(&[true, true], 2, &mut bam_file, Some(Format::Bam)).unwrap();
        }

        let source = AlignmentSource::new(&bam_path);
        let reads_to_keep = vec![false, true];
        let nb_reads_keep = 1;
        let mut buffer = Vec::new();
        
        let result = source.filter_reads_into(&reads_to_keep, nb_reads_keep, &mut buffer, Some(Format::Sam));
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), 10);

        let output = String::from_utf8(buffer).unwrap();
        assert!(output.contains("r02"));
        assert!(output.contains("GCTGACTGAC"));
        assert!(output.contains("JJJJJJJJJJ"));
    }
}
