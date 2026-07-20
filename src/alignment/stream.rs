use std::collections::{BinaryHeap, HashSet};

use anyhow::{Context, Result};
use log::{info, warn};

use noodles::sam::alignment::RecordBuf;
use noodles::sam::Header;
use noodles_util::alignment;
use rand::Rng;

use super::args::Alignment;
use super::io::AlignmentWriter;
use super::model::ScoredRead;
use super::util::extract_name;

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

        let mut survivor_names: HashSet<Vec<u8>> = HashSet::new();

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
        survivor_names: &mut HashSet<Vec<u8>>,
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
        let mut active_reads: BinaryHeap<ScoredRead> = BinaryHeap::with_capacity(target_depth);

        // we also pre alocated a vector here, and will be reused
        let mut survivors: Vec<ScoredRead> = Vec::with_capacity(target_depth);

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
                    for scored_read in &active_reads {
                        // before we write into the result, we need to extract the read name
                        if paired_mode {
                            let qname: Vec<u8> = extract_name(&scored_read.record);

                            survivor_names.insert(qname);
                        }
                        writer
                            .write_record(header, &scored_read.record)
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

            for scored_read in survivors.drain(..) {
                // inclusive end
                if scored_read.end <= start {
                    // the reads survived!

                    // before we write into the result, we need to extract the read name
                    if paired_mode {
                        let qname: Vec<u8> = extract_name(&scored_read.record);

                        survivor_names.insert(qname);
                    }
                    // and we write it into the result
                    writer
                        .write_record(header, &scored_read.record)
                        .context("Failed to write record")?;
                } else {
                    // still active,
                    // put it back in heap
                    active_reads.push(scored_read);
                }
            }

            // note: survivors.drain(..) automatically clears the vector items but keeps the capacity.
            // with target depth usually 50 - 100 maybe, i think this operation is not very expensive.

            if active_reads.len() < target_depth {
                // calculate the end position only once here:
                let end = record
                    .alignment_end()
                    .transpose()?
                    .map(|e| usize::from(e) as i64)
                    .unwrap_or(pos);

                // we convert to Recordbuf here before we want to push it to the heap
                let record = RecordBuf::try_from_alignment_record(header, &record)?;
                let new_read = ScoredRead::new(record, key, end);

                // if the heap has records fewer than N reads, just accept it.
                active_reads.push(new_read);
            } else {
                // if the heap is full, we only consider the new read if it has a better (lower) key
                // than the current worst read in the haep
                if let Some(mut worst) = active_reads.peek_mut() {
                    // does the new record has better key?
                    if key < worst.key {
                        // calculate the distance between these two reads
                        let worst_start = worst.record.alignment_start();
                        let worst_pos = worst_start.map(|p| usize::from(p) as i64 - 1).unwrap_or(0);

                        let distance = pos.saturating_sub(worst_pos);

                        // are we allowed to swap these two reads given a specified swap distance?
                        if distance <= self.swap_distance {
                            // we swap the records
                            // calculate the end position only once here:
                            let end = record
                                .alignment_end()
                                .transpose()?
                                .map(|e| usize::from(e) as i64)
                                .unwrap_or(pos);

                            // we convert to Recordbuf here before we want to swap it to the heap
                            let record = RecordBuf::try_from_alignment_record(header, &record)?;
                            let new_read = ScoredRead::new(record, key, end);
                            // we swap the records
                            *worst = new_read;
                        }
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
        for scored_read in active_reads {
            // before we write into the result, we need to extract the read name
            if paired_mode {
                let qname: Vec<u8> = extract_name(&scored_read.record);

                survivor_names.insert(qname);
            }
            writer
                .write_record(header, &scored_read.record)
                .context("Failed to write record")?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::super::args::SubsamplingStrategy;
    use super::*;
    use noodles_util::alignment::io::Format;
    use std::path::PathBuf;
    use tempfile::NamedTempFile;

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
