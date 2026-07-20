use std::collections::BinaryHeap;

use anyhow::{Context, Result};
use log::{info, warn};

use noodles::sam::Header;
use noodles_util::alignment;
use rand::Rng;

use super::args::Alignment;
use super::io::AlignmentWriter;
use super::model::MappedRead;
use super::NameSet;

/// Whether an active read has fully expired at the current scan position, meaning it can leave
/// the active heap and be written out (its end lies at or before the current start).
fn has_expired(mapped_read: &MappedRead, current_pos: i64) -> bool {
    mapped_read.end <= current_pos
}

/// Whether a new candidate should evict the current worst (highest-key) read in a full heap: the
/// candidate needs a better (lower) priority key, and must lie within the allowed swap distance
/// of the read it would replace.
fn should_swap(
    worst: &MappedRead,
    candidate_key: u64,
    candidate_pos: i64,
    swap_distance: i64,
) -> bool {
    candidate_key < worst.key && candidate_pos.saturating_sub(worst.start) <= swap_distance
}

impl Alignment {
    pub(super) fn run_stream(&self) -> Result<()> {
        info!("Subsampling alignment file (Stream): {:?}", self.aln);

        let mut reader = alignment::io::reader::Builder::default()
            .build_from_path(&self.aln)
            .context("Failed to read alignment file")?;

        let header = reader.read_header()?;

        let (mut rng, mut writer) = self.setup_resources(&header)?;

        // check if reads come from paired end sequencing
        let is_paired = self.check_pair()?;
        let target_depth = self.get_target_depth(is_paired) as usize;

        let mut survivor_names: NameSet =
            NameSet::with_capacity_and_hasher(target_depth, Default::default());

        if is_paired {
            info!("Detected Paired-End data");
        }

        // run sweep line, but only care about first segment record
        // we consider something like nanopore reads are single end, which only have first segments
        self.stream(
            target_depth,
            is_paired,
            &mut survivor_names,
            &header,
            &mut writer,
            &mut rng,
        )?;

        if is_paired {
            // iterate again to recover their mates (or last segement record)
            self.recover_mates(&mut survivor_names, &header, &mut writer)?;
        }
        Ok(())
    }

    // subsampling using sweep line algorithm
    fn stream(
        &self,
        target_depth: usize,
        paired_mode: bool,
        survivor_names: &mut NameSet,
        header: &Header,
        writer: &mut AlignmentWriter,
        rng: &mut rand_pcg::Pcg64,
    ) -> Result<()> {
        // set up reader
        let mut reader = alignment::io::reader::Builder::default()
            .build_from_path(&self.aln)
            .context("Failed to read alignment file")?;
        let _ = reader.read_header()?; // skip header

        // the active set to store reads that currently in the scan
        // using max heap, it keeps the read with the highest key (worst priority) at the top.
        let mut active_reads: BinaryHeap<MappedRead> = BinaryHeap::with_capacity(target_depth);

        // we also pre alocated a vector here, and will be reused
        let mut survivors: Vec<MappedRead> = Vec::with_capacity(target_depth);

        let mut current_tid: Option<usize> = None;
        let mut max_observed_depth: usize = 0;
        // track the last position so that we know whether the input is already sorted or no
        let mut last_pos: i64 = -1;

        // to get chromosome name for warning feature
        let chrom_names: Vec<String> = header
            .reference_sequences()
            .keys()
            .map(|k| k.to_string())
            .collect();

        // a helper closure for giving a warning or report after finished scanning a chromosome
        let depth_report = |tid: Option<usize>, max_depth: usize| {
            if let Some(id) = tid {
                let name = chrom_names.get(id).map(|s| s.as_str()).unwrap_or("*");
                if max_depth < target_depth {
                    if paired_mode {
                        warn!(
                            "Chromosome {} may not have reached target depth (Max: ~{}X [Based on the first segments])",
                            name, max_depth * 2
                        );
                    } else {
                        warn!(
                            "Chromosome {} never reached target depth (Max: {}X)",
                            name, max_depth
                        );
                    }
                } else {
                    info!("Subsampling complete for: {}", name);
                }
            }
        };

        // iterate thriugh every read in the file
        for result in reader.records(header) {
            let record = result.context("Failed to parse BAM record")?;

            // if the read is unmapped, we skip it
            if record.flags()?.is_unmapped() {
                continue;
            }

            // if it is paired, we care only the first segment first
            if paired_mode && record.flags()?.is_last_segment() {
                continue;
            }

            let Some(alignment_start) = record.alignment_start().transpose()? else {
                continue;
            };

            let tid = record.reference_sequence_id(header).transpose()?;

            // noodles use 1-base coordinate.. inclusive start
            let pos = usize::from(alignment_start) as i64 - 1;

            // check if the input is not sorted yet
            let chrom_name = if let Some(id) = tid {
                chrom_names.get(id).map(|n| n.as_str()).unwrap_or("Unknown")
            } else {
                "Unknown"
            };
            if tid == current_tid && pos < last_pos {
                return Err(anyhow::anyhow!(
                    "Input is not sorted! Found read at pos {} after pos {} on {:?}",
                    pos,
                    last_pos,
                    chrom_name
                ));
            }

            if tid != current_tid {
                // we just finished a chromosome. Did it ever reach target coverage?
                depth_report(current_tid, max_observed_depth);

                // we still have survived reads from the previous chromosome in the heap,
                // we need to write them
                if current_tid.is_some() {
                    for mapped_read in &active_reads {
                        // before we write into the result, we need to extract the read name
                        if paired_mode {
                            survivor_names.insert(mapped_read.name.clone());
                        }
                        writer
                            .write_record(header, &mapped_read.record)
                            .context("Failed to write record")?;
                    }
                }
                // Reset for new chromosome
                current_tid = tid;
                max_observed_depth = 0;
                active_reads.clear();
            }

            let start = pos;

            // assign deterministic random priority
            let key: u64 = rng.next_u64();

            // clear and move everything from the heap to the preallocated vector
            // to check whether the reads in the heap already expire
            // the end pos of any read in the heap is lesser or equal than the current start pos.
            // i use drain method: https://doc.rust-lang.org/std/collections/struct.BinaryHeap.html#method.drain
            survivors.extend(active_reads.drain());

            for mapped_read in survivors.drain(..) {
                if has_expired(&mapped_read, start) {
                    // the reads survived!

                    // before we write into the result, we need to extract the read name
                    if paired_mode {
                        survivor_names.insert(mapped_read.name.clone());
                    }
                    // and we write it into the result
                    writer
                        .write_record(header, &mapped_read.record)
                        .context("Failed to write record")?;
                } else {
                    // still active,
                    // put it back in heap
                    active_reads.push(mapped_read);
                }
            }

            // note: survivors.drain(..) automatically clears the vector items but keeps the capacity.
            // with target depth usually 50 - 100 maybe, i think this operation is not very expensive.

            // calculate the end position only once here, regardless of whether the read is accepted:
            let end = record
                .alignment_end()
                .transpose()?
                .map(|e| usize::from(e) as i64)
                .unwrap_or(pos);

            if active_reads.len() < target_depth {
                // the record is already owned (boxed by the reader), so no extra parsing is
                // needed before it enters the heap
                let new_read = MappedRead::new(record, key, start, end);

                // if the heap has records fewer than N reads, just accept it.
                active_reads.push(new_read);
            } else {
                // if the heap is full, we only consider the new read if it has a better (lower) key
                // than the current worst read in the heap, and lies within the swap distance
                if let Some(mut worst) = active_reads.peek_mut() {
                    if should_swap(&worst, key, pos, self.swap_distance) {
                        let new_read = MappedRead::new(record, key, start, end);
                        *worst = new_read;
                    }
                }
            }

            // update the last position
            last_pos = pos;

            // track stats for Warning
            max_observed_depth = max_observed_depth.max(active_reads.len());
        }

        // report a warning or info forn the last chromosome
        depth_report(current_tid, max_observed_depth);

        // any remaining read in the heap survive
        // write into the subsampled result.
        for mapped_read in active_reads {
            // before we write into the result, we need to extract the read name
            if paired_mode {
                survivor_names.insert(mapped_read.name.clone());
            }
            writer
                .write_record(header, &mapped_read.record)
                .context("Failed to write record")?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::super::args::SubsamplingStrategy;
    use super::*;
    use noodles::sam::alignment::{Record, RecordBuf};
    use noodles_util::alignment::io::Format;
    use std::path::PathBuf;
    use tempfile::NamedTempFile;

    fn mapped_read(key: u64, start: i64, end: i64) -> MappedRead {
        let record: Box<dyn Record> = Box::new(RecordBuf::default());
        MappedRead::new(record, key, start, end)
    }

    #[test]
    fn test_has_expired_when_end_before_current_pos() {
        let read = mapped_read(1, 0, 100);
        assert!(has_expired(&read, 150));
    }

    #[test]
    fn test_has_expired_when_end_equals_current_pos() {
        // end is inclusive, so a read ending exactly at the scan position has expired
        let read = mapped_read(1, 0, 100);
        assert!(has_expired(&read, 100));
    }

    #[test]
    fn test_not_expired_when_end_after_current_pos() {
        let read = mapped_read(1, 0, 100);
        assert!(!has_expired(&read, 99));
    }

    #[test]
    fn test_should_swap_with_better_key_within_distance() {
        let worst = mapped_read(50, 10, 100);
        assert!(should_swap(&worst, 5, 12, 5));
    }

    #[test]
    fn test_should_not_swap_with_worse_key() {
        let worst = mapped_read(5, 10, 100);
        assert!(!should_swap(&worst, 50, 12, 5));
    }

    #[test]
    fn test_should_not_swap_beyond_swap_distance() {
        let worst = mapped_read(50, 10, 100);
        assert!(!should_swap(&worst, 5, 20, 5));
    }

    #[test]
    fn test_should_swap_at_exact_swap_distance_boundary() {
        let worst = mapped_read(50, 10, 100);
        assert!(should_swap(&worst, 5, 15, 5));
    }

    #[test]
    fn test_output_never_exceeds_target_depth() {
        let input_path = PathBuf::from("tests/cases/test.bam");
        let target_depth = 3;

        let output_file = NamedTempFile::new().unwrap();
        let output_path = output_file.path().to_path_buf();

        let align = Alignment {
            aln: input_path,
            output: Some(output_path.clone()),
            output_format: Some(Format::Bam),
            coverage: target_depth,
            seed: Some(2109),
            strategy: SubsamplingStrategy::Stream,
            swap_distance: 5,
            step_size: 100, // it is not used
            batch_size: 10_000,
        };

        // run the sweep line algorithm
        align.run_stream().expect("Subsampling failed");

        // verify the resulted depth
        let mut reader = alignment::io::reader::Builder::default()
            .build_from_path(&output_path)
            .unwrap();

        let header = reader.read_header().unwrap();

        let mut records: Vec<RecordBuf> = reader
            .records(&header)
            .map(|r| RecordBuf::try_from_alignment_record(&header, &r.unwrap()).unwrap())
            .collect();

        // sort by chromosome and then by position
        records.sort_by(|a, b| {
            let tid_a = a.reference_sequence_id().unwrap_or(usize::MAX);
            let tid_b = b.reference_sequence_id().unwrap_or(usize::MAX);
            let pos_a = a.alignment_start().map(usize::from).unwrap_or(0);
            let pos_b = b.alignment_start().map(usize::from).unwrap_or(0);

            match tid_a.cmp(&tid_b) {
                std::cmp::Ordering::Equal => pos_a.cmp(&pos_b),
                other => other,
            }
        });

        let mut active_ends: Vec<usize> = Vec::new();
        let mut max_obs_depth = 0;
        let mut current_tid = None;

        for record in records {
            let tid = record.reference_sequence_id();
            // Convert to 0-based inclusive start for comparison logic
            let start = usize::from(record.alignment_start().unwrap()) - 1;
            // Convert to 0-based exclusive end (which equals 1-based inclusive end value)
            let end = usize::from(record.alignment_end().unwrap());

            if tid != current_tid {
                active_ends.clear();
                current_tid = tid;
            }

            // retain reads that end AFTER the current start position
            active_ends.retain(|&e| e > start);

            active_ends.push(end);

            let current_depth = active_ends.len();
            if current_depth > max_obs_depth {
                max_obs_depth = current_depth;
            }
        }

        assert!(
            max_obs_depth <= target_depth as usize,
            "The resulted depth {} exceeds the target depth {}",
            max_obs_depth,
            target_depth
        );
    }

    #[test]
    fn test_subsampling_statistics() {
        let input_path = PathBuf::from("tests/cases/test.bam");

        let target_depth = 3;

        let output_temp = NamedTempFile::new().unwrap();

        let align = Alignment {
            aln: input_path,
            output: Some(output_temp.path().to_path_buf()),
            output_format: Some(Format::Bam),
            coverage: target_depth,
            seed: Some(2109),
            strategy: SubsamplingStrategy::Stream,
            swap_distance: 5,
            step_size: 100, // not used
            batch_size: 10_000,
        };
        align.run_stream().expect("Subsampling failed");

        let mut reader = alignment::io::reader::Builder::default()
            .build_from_path(output_temp.path())
            .unwrap();
        // read the header, get the length of chromosome
        let header = reader.read_header().unwrap();

        let chrom_tid = header
            .reference_sequences()
            .get_index_of(b"plasmid_2".as_slice())
            .expect("Chromosome not found");
        let chrom_length = header
            .reference_sequences()
            .get_index(chrom_tid)
            .unwrap()
            .1
            .length()
            .get();

        // using preallocated vector to record depth for each position
        // populate the vector with 0
        let mut depth = vec![0u32; chrom_length];

        for r in reader.records(&header) {
            let record = r.unwrap();
            let r_tid = record.reference_sequence_id(&header).unwrap();

            // only check plasmid_2 contig. It has:
            // min 2, median 7, mean 7.5761410788382, max 11.
            if r_tid.unwrap() != chrom_tid {
                continue;
            }

            // convert to 0 based for vector indexing
            let start = usize::from(record.alignment_start().unwrap().unwrap()) - (1_usize);
            let end = usize::from(record.alignment_end().unwrap().unwrap());
            // update the posisition depth that were covered by the record
            for pos in start..end {
                if (pos) < depth.len() {
                    depth[pos] += 1; // update the depth
                } // safety
            }
        }

        if depth.is_empty() {
            panic!("Resulting BAM is empty! Subsampling failed completely.");
        }

        let sum_depth: u64 = depth.iter().map(|&d| d as u64).sum();
        let mean = sum_depth as f64 / chrom_length as f64;

        // sort to find median
        depth.sort_unstable();
        let mid = chrom_length / 2;
        let median = depth[mid];

        // check if the median is exactly the target depth or at least close to
        assert!(
            median >= target_depth - 1 && median <= target_depth,
            "Median depth {} is too far from target {}",
            median,
            target_depth
        );

        // check if the mean is exactly the target depth or at least close to
        assert!(
            (mean - target_depth as f64).abs() < 1.0,
            "Mean depth {:.2} deviated too much from target {}",
            mean,
            target_depth
        );
    }
}
