use std::cmp::Reverse;
use std::collections::BinaryHeap;

use anyhow::{anyhow, Context, Result};
use log::{info, warn};
use rand::prelude::*;

use noodles::core::{Position, Region};
use noodles::sam::alignment::Record;
use noodles::sam::Header;
use noodles_util::alignment::io::indexed_reader;

use super::args::Alignment;
use super::io::AlignmentWriter;
use super::util::{alignment_start, extract_name, shuffle_grouped_by_position};
use super::NameSet;

impl Alignment {
    pub(super) fn run_fetch(&self) -> Result<()> {
        info!("Subsampling alignment file (Fetch): {:?}", self.aln);

        let mut reader = indexed_reader::Builder::default()
            .build_from_path(&self.aln)
            .context("Failed to read alignment file")?;

        let header = reader.read_header()?;

        let (mut rng, mut writer) = self.setup_resources(&header)?;

        let is_paired = self.check_pair()?;

        let target_depth = self.get_target_depth(is_paired); // no memory allocation, do not need to convert into usize
        let mut survivor_names: NameSet =
            NameSet::with_capacity_and_hasher(target_depth as usize, Default::default());

        if is_paired {
            info!("Detected Paired-End data");
        }

        // run fetching but only care about first segment
        // we consider something like nanopore reads are single end, which only have first segments
        self.fetching(
            target_depth,
            is_paired,
            &mut survivor_names,
            &header,
            &mut writer,
            &mut rng,
        )?;

        if is_paired {
            // recover their mates
            self.recover_mates(&mut survivor_names, &header, &mut writer)?;
        }

        Ok(())
    }

    fn fetching(
        &self,
        target_depth: u32,
        paired_mode: bool,
        survivor_names: &mut NameSet,
        header: &Header,
        writer: &mut AlignmentWriter,
        rng: &mut rand_pcg::Pcg64,
    ) -> Result<()> {
        let mut reader = indexed_reader::Builder::default()
            .build_from_path(&self.aln)
            .context("Failed to read alignment file")?;

        let _ = reader.read_header()?;

        // to hold batch of reads inside batch size query. Kept as the boxed record the reader
        // already owns rather than eagerly parsed into a `RecordBuf` -- most cached reads are
        // never selected, so paying for a full parse of every field of every read in the batch
        // would be wasted work.
        let mut batch_cache: Vec<Box<dyn Record>> = Vec::new();

        // how big of a chunk to grab at once
        let batch_size = self.batch_size as usize;

        // iterate through directly over reference sequences
        for (tid, (chrom_name, _)) in header.reference_sequences().iter().enumerate() {
            let chrom_name_str = String::from_utf8_lossy(chrom_name.as_ref());
            info!("Subsampling chromosome: {chrom_name}");

            // get chromosome length
            let chrom_len = header
                .reference_sequences()
                .get_index(tid)
                .map(|(_, rs)| rs.length().get() as u64)
                .context(format!(
                    "Failed to get chromosome length for chromosome {chrom_name}"
                ))?;

            // get first record position
            let first_pos = {
                let full_region = Region::new(chrom_name.clone(), ..);
                let mut query = reader.query(header, &full_region)?;
                // Peek the first record to find start position
                match query.next() {
                    Some(Ok(r)) => r
                        .alignment_start()
                        .transpose()?
                        .map(usize::from)
                        .unwrap_or(1),
                    _ => {
                        warn!("Chromosome {} has no records", chrom_name_str);
                        continue;
                    }
                }
            };

            let mut next_pos = first_pos;
            let start_pos = next_pos;

            let mut n_reads_needed = target_depth;
            let mut current_reads: NameSet = NameSet::default();
            let mut heap: BinaryHeap<Reverse<(usize, Vec<u8>)>> = BinaryHeap::new();
            let mut regions_below_coverage = false;

            // clear the cache for new chromosome
            batch_cache.clear();

            // the region that is currently covered
            let mut cache_start = 0;
            let mut cache_end = 0;

            loop {
                // noodles Position::new takes 1-based index

                // check if next_pos is outside the current cache
                if next_pos < cache_start || next_pos >= cache_end {
                    // refill the cache
                    cache_start = next_pos;
                    cache_end = next_pos + batch_size;

                    // get a whole chunk of reads in this region
                    let start = Position::new(next_pos).unwrap_or(Position::MIN);
                    // ensure the end of region doesn't exceed the length of chrom
                    let safe_end = (cache_end).min(chrom_len as usize);
                    let end = Position::new(safe_end).unwrap_or(start);
                    let region = Region::new(chrom_name_str.as_ref(), start..=end);

                    // clear it for new batch
                    batch_cache.clear();
                    // query reads that intersect the region
                    let query = reader.query(header, &region)?;

                    // put them in the cache. The reader already hands back owned records, so no
                    // conversion is needed here -- only the ones actually selected below get
                    // written out.
                    for result in query {
                        let record = result?;
                        batch_cache.push(record);
                    }
                }

                // get reads that overlap with next_pos (shuffle point) from batch (already in the memory owned)
                // similar as `reader.query(next_pos)`

                // reads that overlaps at a certain position
                let mut candidates: Vec<&Box<dyn Record>> = Vec::new();
                // ignore records that start after next_pos, https://doc.rust-lang.org/std/primitive.slice.html#method.partition_point
                // i assume the input is sorted, because it requires the index
                let partition_idx = batch_cache.partition_point(|r| {
                    let start = alignment_start(r).unwrap_or(Position::MIN);
                    usize::from(start) <= next_pos
                });

                for record in &batch_cache[..partition_idx] {
                    if paired_mode && record.flags()?.is_last_segment() {
                        continue;
                    }

                    let end =
                        usize::from(record.alignment_end().transpose()?.unwrap_or(Position::MIN));

                    // just check if next_pos "inside" the reads
                    if end >= next_pos {
                        candidates.push(record);
                    }
                }

                if next_pos == start_pos {
                    candidates.shuffle(rng);
                } else {
                    // the batch is already sorted by position (it came straight from an indexed
                    // query), so we only need to group-shuffle, not re-sort it first
                    shuffle_grouped_by_position(&mut candidates, rng);
                }

                let mut num_output = 0;
                let mut record_iter = candidates.iter().rev();

                while num_output < n_reads_needed {
                    let record = match record_iter.next() {
                        Some(r) => *r,
                        None => break,
                    };

                    let qname: Vec<u8> = extract_name(record);

                    if current_reads.contains(&qname) {
                        continue;
                    }

                    let end = record
                        .alignment_end()
                        .transpose()?
                        .map(usize::from)
                        .unwrap_or(next_pos);

                    current_reads.insert(qname.clone());
                    heap.push(Reverse((end, qname.clone())));

                    // save record names, if it is paired end data
                    if paired_mode {
                        survivor_names.insert(qname);
                    }

                    // write the record
                    writer
                        .write_record(header, record)
                        .context("Failed to write record")?;
                    num_output += 1;
                }

                n_reads_needed = n_reads_needed.saturating_sub(num_output);

                if n_reads_needed > 0 {
                    // increment next_pos by step_size or the minimum end position of the reads in
                    // the region, whichever is smaller
                    let min_end = heap.peek().map(|Reverse((end, _))| *end);
                    let jump = match min_end {
                        Some(min) => (self.step_size as usize).min(min.saturating_sub(next_pos)),
                        None => self.step_size as usize,
                    };
                    next_pos += jump;
                    regions_below_coverage = true;
                }

                // remove smallest end position from heap and update next_pos and n_reads_needed
                // allowing for the fact that there may be multiple reads with the same end position
                while let Some(Reverse((end, qname))) = heap.pop() {
                    // because noodles using 1-based index inclusive
                    next_pos = end + 1;

                    if current_reads.contains(&qname) {
                        current_reads.remove(&qname);
                        n_reads_needed += 1;
                    } else {
                        return Err(anyhow!("A read in the heap was not found in the current reads set. This should not happen, please raise an issue."));
                    }

                    if heap.is_empty() {
                        break;
                    }
                    if let Some(Reverse((next_end, _))) = heap.peek() {
                        if *next_end != end {
                            break;
                        }
                    } else {
                        break;
                    }
                }
                if next_pos >= chrom_len as usize {
                    break;
                }
            }
            if regions_below_coverage {
                if paired_mode {
                    warn!(
                        "Chromosome {chrom_name} has regions that may have less than the requested coverage [Based on the first segments]");
                } else {
                    warn!(
                        "Chromosome {chrom_name} has regions with less than the requested coverage"
                    );
                }
            }
        }
        Ok(())
    }
}
