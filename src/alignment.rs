use std::borrow::Cow;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::result::Result::Ok;

use anyhow::{anyhow, Context, Result};
use clap::{Parser, ValueEnum};
use log::{info, warn};
use rand::prelude::SliceRandom;
use rand::{random, Rng, SeedableRng};

use noodles::core::{Position, Region};
use noodles::sam::alignment::{Record, RecordBuf};
use noodles::sam::header::record::value::{
    map::{program::tag, Program},
    Map,
};

use noodles::sam::Header;
use noodles_util::alignment::{
    self,
    io::{indexed_reader, Format},
};

// a generic writer that can handle File or Stdout
type AlignmentWriter = alignment::io::Writer<Box<dyn Write>>;

use crate::cli::check_path_exists;
use crate::Runner;

const RASUSA: &str = "rasusa";

// we implementes traits manually to avoid comparing the full record
// primary comparison is based using the random priority 'key'
// if there is tie (two read have the same key), then we use 'qname' to ensure the deterministic sorting still
#[derive(Debug)]
struct ScoredRead {
    /// Alignment record (SAM/BAM/CRAM)
    record: RecordBuf,

    /// The deterministic key assigned to this record
    key: u64,

    /// Alignment end for the record (1-based inclusive)
    end: i64,
}

impl ScoredRead {
    fn new(record: RecordBuf, key: u64, end: i64) -> Self {
        Self { record, key, end }
    }
}

impl PartialEq for ScoredRead {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key && self.record.name() == other.record.name()
    }
}

impl Eq for ScoredRead {}

impl PartialOrd for ScoredRead {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

// if the random key is equal (very rare),
// then compare it based the qname record (which i think will be unique)
impl Ord for ScoredRead {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.key.cmp(&other.key) {
            std::cmp::Ordering::Equal => self.record.name().cmp(&other.record.name()),
            std::cmp::Ordering::Greater => std::cmp::Ordering::Greater,
            std::cmp::Ordering::Less => std::cmp::Ordering::Less,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, ValueEnum)]
pub enum SubsamplingStrategy {
    /// A linear scan approach using sweep line algorithm with random priority. Requires sorted alignment input.
    Stream,

    /// A fetching approach to randomly subsample reads given read overlap position. Requires indexed input (.bai).
    Fetch,
}

#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct Alignment {
    /// Path to the input alignment file (SAM/BAM/CRAM) to subsample
    ///
    /// Note: An index (.bai) is required when using '--strategy fetch'.
    #[arg(value_parser = check_path_exists, name = "FILE")]
    pub aln: PathBuf,

    /// Path to the output subsampled alignment file. Defaults to stdout (same format as input)
    ///
    /// The output is not guaranteed to be sorted. We recommend piping the output to `samtools sort`
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,

    /// Output format. Rasusa will attempt to infer the format from the output file extension if not provided
    #[arg(short='O', long, value_name = "FMT", value_parser = infer_format_from_char)]
    pub output_type: Option<Format>,

    /// The desired depth of coverage to subsample the alignment to
    #[arg(short, long, value_name = "INT", value_parser = clap::value_parser!(u32).range(1..))]
    pub coverage: u32,

    /// Random seed to use.
    #[arg(short, long, value_name = "INT")]
    pub seed: Option<u64>,

    // Strategy selection
    /// Subsampling strategy
    #[arg(long, value_enum, default_value_t = SubsamplingStrategy::Stream)]
    pub strategy: SubsamplingStrategy,

    // Algorithm specific arguments
    /// [Stream] A maximum distance (bp) allowed between start position of new read and the worst read
    /// in the heap to consider them to be 'swappable'.
    ///
    /// Larger values allow swapping reads over greater distances, but may cause local undersampling.
    /// A value of `0` means only allows swap between reads that have the same start position.
    #[arg(long, default_value_t = 5, value_name = "INT", value_parser = clap::value_parser!(i64).range(0..))]
    pub swap_distance: i64,

    /// [Fetch] When a region has less than the desired coverage, the step size to move along the chromosome
    /// to find more reads.
    ///
    /// The lowest of the step and the minimum end coordinate of the reads in the region will be used.
    /// This parameter can have a significant impact on the runtime of the subsampling process.
    #[arg(long, default_value_t = 100, value_name = "INT", value_parser = clap::value_parser!(i64).range(1..))]
    pub step_size: i64,

    /// [Fetch] The size of the genomic window (bp) to cache into memory at once.
    ///
    /// Larger values reduce disk seeking, but at the cost of high memory usage.
    /// The minimum value is 1,000 bp to avoid small region queries.
    #[arg(long, default_value_t = 10_000, value_name = "INT", value_parser = clap::value_parser!(u64).range(1000..))]
    pub batch_size: u64,
}

impl Runner for Alignment {
    fn run(&mut self) -> Result<()> {
        match self.strategy {
            SubsamplingStrategy::Stream => self.run_stream(),
            SubsamplingStrategy::Fetch => self.run_fetch(),
        }
    }
}

impl Alignment {

    fn setup_resources(&self, input_header: &Header) -> Result<(rand_pcg::Pcg64, AlignmentWriter)> {
        // set up random number generator
        let rng = match self.seed {
            Some(s) => rand_pcg::Pcg64::seed_from_u64(s),
            None => {
                let seed = random();
                info!("Using seed: {seed}");
                rand_pcg::Pcg64::seed_from_u64(seed)
            }
        };

        let mut header = input_header.clone();

        // add rasusa program command line to header
        let (pg_id, pg_map) = self.program_entry(&header);

        // set up the header and writer
        header.programs_mut().as_mut().insert(pg_id.into(), pg_map);

        let input_fmt = match infer_format_from_path(&self.aln) {
            Some(fmt) => fmt,
            None => {
                return Err(anyhow::anyhow!(
                    "Output file format not recognized. Please use .sam, .bam, or .cram extensions"
                ));
            }
        };

        let output_fmt = match &self.output_type {
            Some(fmt) => *fmt,
            None => match &self.output {
                None => input_fmt,
                Some(path) => match infer_format_from_path(path) {
                    Some(fmt) => fmt,
                    None => {
                        return Err(anyhow::anyhow!(
                            "Output file format not recognized. Please use .sam, .bam, or .cram extensions"
                        ));
                    }
                },
            },
        };

        // use Box<dyn Write> to make File and Stdout compatible
        let sink: Box<dyn Write> = match &self.output {
            Some(path) => {
                let path = path.as_path();
                info!("Writing subsampled alignment to: {:?}", path);
                let file = File::create(path).context("Failed to create output alignment file")?;
                Box::new(file)
            }
            None => {
                info!("Writing subsampled alignment to stdout");
                Box::new(io::stdout().lock())
            }
        };

        let mut writer = alignment::io::writer::Builder::default()
            .set_format(output_fmt)
            .build_from_writer(sink)?;

        writer.write_header(&header)?;

        Ok((rng, writer))
    }

    // a function which infers data type (paired and or single end) by looking at the first 10 records
    fn check_pair(&self) -> Result<bool> {
        // set up reader
        let mut reader = alignment::io::reader::Builder::default()
            .build_from_path(&self.aln)
            .context("Failed to read alignment file")?;

        let header = reader.read_header()?;

        // read records only to check what type of the data
        let mut check_n: u8 = 1;
        for result in reader.records(&header) {
            let record = result.context("Failed to parse BAM record")?;
            if record.flags()?.is_segmented() {
                return Ok(true);
            }
            if check_n > 10 {
                break;
            }
            check_n += 1;
        }
        Ok(false)
    }

    // a helper function to calculate target depth based on their input
    fn get_target_depth(&self, is_paired: bool) -> u32 {
        if is_paired {
            // divide the original target by 2 because we will add the mates back later
            (self.coverage / 2).max(1) // in case we hit 0
        } else {
            self.coverage
        }
    }

    // a helper function to scan the file linearly to find mates, only relevent on paired end illumina data
    fn recover_mates(
        &self,
        survivor_names: &mut HashSet<Vec<u8>>,
        header: &Header,
        writer: &mut AlignmentWriter,
    ) -> Result<()> {
        info!("Recovering mates (last segment records)");
        let mut reader = alignment::io::reader::Builder::default()
            .build_from_path(&self.aln)
            .context("Failed to read alignment file")?;

        let _ = reader.read_header()?; // skip header

        for result in reader.records(header) {
            let record = result.context("Failed to parse BAM record")?;

            // only care about the last segment (read 2)
            if !record.flags()?.is_last_segment() {
                continue;
            }

            let qname: Vec<u8> = extract_name(&record);

            if survivor_names.contains(&qname) {
                writer
                    .write_record(header, &record)
                    .context("Failed to write mate records")?;
            }
        }
        Ok(())
    }

    fn run_stream(&self) -> Result<()> {
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
            let key: u64 = rng.gen();

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

    fn run_fetch(&self) -> Result<()> {
        info!("Subsampling alignment file (Fetch): {:?}", self.aln);

        let mut reader = indexed_reader::Builder::default()
            .build_from_path(&self.aln)
            .context("Failed to read alignment file")?;

        let header = reader.read_header()?;

        let (mut rng, mut writer) = self.setup_resources(&header)?;

        let is_paired = self.check_pair()?;

        let target_depth = self.get_target_depth(is_paired); // no memory allocation, do not need to convert into usize
        let mut survivor_names: HashSet<Vec<u8>> = HashSet::new();

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
        survivor_names: &mut HashSet<Vec<u8>>,
        header: &Header,
        writer: &mut AlignmentWriter,
        rng: &mut rand_pcg::Pcg64,
    ) -> Result<()> {
        let mut reader = indexed_reader::Builder::default()
            .build_from_path(&self.aln)
            .context("Failed to read alignment file")?;

        let _ = reader.read_header()?;

        // to hold batch of reads inside batch size query
        let mut batch_cache: Vec<RecordBuf> = Vec::new();

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
            let mut current_reads = HashSet::new();
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

                    // put them in the cache
                    for result in query {
                        let record = result?;
                        // convert to owned Record, because we want shuffle these reads
                        let recordbuf = RecordBuf::try_from_alignment_record(header, &record)?;
                        batch_cache.push(recordbuf);
                    }
                }

                // get reads that overlap with next_pos (shuffle point) from batch (already in the memory owned)
                // similar as `reader.query(next_pos)`

                // reads that overlaps at a certain position
                let mut candidates: Vec<&RecordBuf> = Vec::new();
                // ignore recordbuf that start after next_pos, https://doc.rust-lang.org/std/primitive.slice.html#method.partition_point
                // i assume the input is sorted, because it requires the index
                let partition_idx = batch_cache.partition_point(|r| {
                    usize::from(r.alignment_start().unwrap_or(Position::MIN)) <= next_pos
                });

                for record in &batch_cache[..partition_idx] {
                    if paired_mode && record.flags().is_last_segment() {
                        continue;
                    }

                    let end = usize::from(record.alignment_end().unwrap_or(Position::MIN));

                    // just check if next_pos "inside" the reads
                    if end >= next_pos {
                        candidates.push(record);
                    }
                }

                if next_pos == start_pos {
                    candidates.shuffle(rng);
                } else {
                    shuffle_records_by_position(&mut candidates, rng);
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

                    let end = record.alignment_end().map(usize::from).unwrap_or(next_pos);

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

    /// Generates a rasusa program entry from a SAM header
    fn program_entry(&self, header: &Header) -> (String, Map<Program>) {
        let (program_id, previous_pgid) = make_program_id_unique(header, RASUSA);

        // Creates a SAM header record map value
        let mut record = Map::<Program>::builder();

        record = record.insert(tag::NAME, RASUSA);
        record = record.insert(tag::VERSION, env!("CARGO_PKG_VERSION"));

        let cl = std::env::args().collect::<Vec<String>>().join(" ");
        record = record.insert(tag::COMMAND_LINE, cl);

        // Link to previous program
        if let Some(pp) = previous_pgid {
            record = record.insert(tag::PREVIOUS_PROGRAM_ID, pp);
        };

        let program = record.build().expect("Failed to build program record");

        (program_id.into_owned(), program)
    }
}

/// Sorts records by position and shuffles those with the same position
fn shuffle_records_by_position(records: &mut [&RecordBuf], rng: &mut impl Rng) {
    // First sort by position
    records.sort_by_key(|record| record.alignment_start());

    // Then shuffle groups with the same position
    let mut start = 0;
    while start < records.len() {
        let pos = records[start].alignment_start();
        let mut end = start + 1;

        // Find the end of the group with the same position
        while end < records.len() && records[end].alignment_start() == pos {
            end += 1;
        }

        // Shuffle this group if it has more than one element
        if end - start > 1 {
            records[start..end].shuffle(rng);
        }

        start = end;
    }
}

/// Makes a program ID unique by looking for existing program records with the same ID and adding
/// a suffix to the ID if necessary. Also returns the program ID of the last program in the header
fn make_program_id_unique<'a>(
    header: &Header,
    program_id: &'a str,
) -> (Cow<'a, str>, Option<String>) {
    let programs = header.programs().as_ref();

    // noodles uses an IndexMap, so .last() will guaranteed to be the most recent entry
    let last_pg_id = programs.keys().last().map(|pp| pp.to_string());

    // count occurance
    let occurrences_of_id = programs
        .keys()
        .filter(|pp| {
            let id = pp.to_string();

            // split "rasusa.1" ->"rasusa"
            let id_before_last_dot = id.rfind('.').map(|i| &id[..i]).unwrap_or(&id);

            id_before_last_dot == program_id
        })
        .count();

    if occurrences_of_id == 0 {
        (Cow::Borrowed(program_id), last_pg_id)
    } else {
        let new_id = format!("{program_id}.{occurrences_of_id}");
        (Cow::Owned(new_id), last_pg_id)
    }
}

// a generic function to extract a read name as bytes
fn extract_name<R: Record>(record: &R) -> Vec<u8> {
    record
        .name()
        .map(|name| {
            let b: &[u8] = name.as_ref();
            b.to_vec()
        })
        .unwrap_or_default()
}

pub fn infer_format_from_path(path: &Path) -> Option<Format> {
    match path.extension().and_then(|ext| ext.to_str()) {
        Some("sam") => Some(Format::Sam),
        Some("bam") => Some(Format::Bam),
        Some("cram") => Some(Format::Cram),
        _ => None,
    }
}

// a function which infers a format from a single character (case insensitive) this will be used in the CLI, so should return a Result
pub fn infer_format_from_char(c: &str) -> Result<Format, String> {
    match c.to_ascii_lowercase().as_str() {
        "s" | "sam" => Ok(Format::Sam),
        "b" | "bam" => Ok(Format::Bam),
        "c" | "cram" => Ok(Format::Cram),
        _ => Err(String::from("Invalid output format. Please use 's', 'b', or 'c' to specify SAM, BAM, or CRAM format, respectively")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_cmd::Command;
    use noodles::core::Position;
    use noodles::sam::alignment::{Record, RecordBuf};
    use rand::prelude::StdRng;
    use tempfile::NamedTempFile;

    const SUB: &str = "aln";

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
    fn test_shuffle_records_by_position_empty() {
        let mut rng = StdRng::seed_from_u64(1234);
        let mut empty_records: Vec<&RecordBuf> = vec![];

        shuffle_records_by_position(&mut empty_records, &mut rng);

        assert_eq!(empty_records.len(), 0);
    }

    #[test]
    fn test_shuffle_records_by_position_single_record() {
        let mut rng = StdRng::seed_from_u64(1234);

        // Create a single BAM record
        let mut record = RecordBuf::default();
        *record.alignment_start_mut() = Position::new(100);
        let mut records: Vec<&RecordBuf> = vec![&record];

        shuffle_records_by_position(&mut records, &mut rng);

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].alignment_start(), Position::new(100));
    }

    #[test]
    fn test_shuffle_records_by_position_maintains_sort_order() {
        let mut rng = StdRng::seed_from_u64(1234);

        // Create records with different positions
        let mut record1 = RecordBuf::default();
        *record1.alignment_start_mut() = Position::new(300);
        let mut record2 = RecordBuf::default();
        *record2.alignment_start_mut() = Position::new(100);
        let mut record3 = RecordBuf::default();
        *record3.alignment_start_mut() = Position::new(200);

        let mut records = vec![&record1, &record2, &record3];

        shuffle_records_by_position(&mut records, &mut rng);

        // Should be sorted by position
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].alignment_start(), Position::new(100));
        assert_eq!(records[1].alignment_start(), Position::new(200));
        assert_eq!(records[2].alignment_start(), Position::new(300));
    }

    #[test]
    fn test_shuffle_records_by_position_shuffles_same_position() {
        let mut rng = StdRng::seed_from_u64(1234);

        // Create multiple records with the same position but different data
        let mut record1 = RecordBuf::default();
        *record1.alignment_start_mut() = Position::new(100);
        *record1.name_mut() = Some("read1".parse().unwrap());

        let mut record2 = RecordBuf::default();
        *record2.alignment_start_mut() = Position::new(100);
        *record2.name_mut() = Some("read2".parse().unwrap());

        let mut record3 = RecordBuf::default();
        *record3.alignment_start_mut() = Position::new(100);
        *record3.name_mut() = Some("read3".parse().unwrap());

        let records = vec![&record1, &record2, &record3];
        let original_order: Vec<String> = records
            .iter()
            .map(|r| r.name().unwrap().to_string())
            .collect();

        // Run shuffle multiple times to verify randomness
        let mut different_orders = 0;
        for _ in 0..10 {
            let mut test_records = records.clone();
            shuffle_records_by_position(&mut test_records, &mut rng);

            // All should still be at position 100
            assert!(records
                .iter()
                .all(|r| r.alignment_start() == Position::new(100)));

            // Check if order changed
            let new_order: Vec<String> = test_records
                .iter()
                .map(|r| r.name().unwrap().to_string())
                .collect();
            if new_order != original_order {
                different_orders += 1;
            }
        }

        // Should have at least some different orders due to shuffling
        assert!(
            different_orders > 0,
            "Records with same position should be shuffled"
        );
    }

    #[test]
    fn test_shuffle_records_by_position_mixed_positions() {
        let mut rng = StdRng::seed_from_u64(1234);

        // Create records with mixed positions - some same, some different
        let mut record1 = RecordBuf::default();
        *record1.alignment_start_mut() = Position::new(100);
        *record1.name_mut() = Some("read1_pos100".parse().unwrap());

        let mut record2 = RecordBuf::default();
        *record2.alignment_start_mut() = Position::new(200);
        *record2.name_mut() = Some("read2_pos200".parse().unwrap());

        let mut record3 = RecordBuf::default();
        *record3.alignment_start_mut() = Position::new(100);
        *record3.name_mut() = Some("read3_pos100".parse().unwrap());

        let mut record4 = RecordBuf::default();
        *record4.alignment_start_mut() = Position::new(150);
        *record4.name_mut() = Some("read4_pos150".parse().unwrap());

        let mut record5 = RecordBuf::default();
        *record5.alignment_start_mut() = Position::new(100);
        *record5.name_mut() = Some("read5_pos100".parse().unwrap());

        let mut records = vec![&record1, &record2, &record3, &record4, &record5];

        shuffle_records_by_position(&mut records, &mut rng);

        // Should be sorted by position
        let positions: Vec<Position> = records
            .iter()
            .map(|r| r.alignment_start().unwrap())
            .collect();
        assert_eq!(
            positions,
            vec![
                Position::new(100).unwrap(),
                Position::new(100).unwrap(),
                Position::new(100).unwrap(),
                Position::new(150).unwrap(),
                Position::new(200).unwrap()
            ]
        );

        // Records at position 100 should potentially be in different order
        let pos100_names: Vec<Vec<u8>> = records
            .iter()
            .filter(|r| r.alignment_start() == Position::new(100))
            .map(|r| r.name().unwrap().to_vec())
            .collect();
        assert_eq!(pos100_names.len(), 3);

        // All names should be present
        let name_strings: Vec<String> = pos100_names
            .iter()
            .map(|name| String::from_utf8_lossy(name).to_string())
            .collect();
        assert!(name_strings.contains(&"read1_pos100".to_string()));
        assert!(name_strings.contains(&"read3_pos100".to_string()));
        assert!(name_strings.contains(&"read5_pos100".to_string()));
    }

    #[test]
    fn no_coverage_given_raises_error() {
        let infile = "tests/cases/test.bam";
        let passed_args = vec![SUB, infile];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn zero_coverage_raises_error_stream() {
        let infile = "tests/cases/test.bam";
        let passed_args = vec![SUB, infile, "-c", "0"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn zero_coverage_raises_error_fetch() {
        let infile = "tests/cases/test.bam";
        let passed_args = vec![
            SUB,
            infile,
            "-c",
            "0",
            "--strategy",
            "fetch",
            "--step-size",
            "5000",
        ];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn bam_with_regions_of_zero_coverage_doesnt_endless_loop_stream() {
        let infile = "tests/cases/test.bam";
        let passed_args = vec![SUB, infile, "-c", "1"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().success();
    }

    #[test]
    fn bam_with_regions_of_zero_coverage_doesnt_endless_loop_fetch() {
        let infile = "tests/cases/test.bam";
        let passed_args = vec![SUB, infile, "-c", "1", "--strategy", "fetch"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().success();
    }

    #[test]
    fn excess_coverage_doesnt_endless_loop_stream() {
        let infile = "tests/cases/test.bam";
        let passed_args = vec![SUB, infile, "-c", "10000"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().success();
    }

    #[test]
    fn excess_coverage_doesnt_endless_loop_fetch() {
        let infile = "tests/cases/test.bam";
        let passed_args = vec![SUB, infile, "-c", "10000", "--strategy", "fetch"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().success();
    }

    // because it doesn't require index anymore
    #[test]
    fn bam_no_index_is_ok_stream() {
        let infile = "tests/cases/no_index.bam";
        let passed_args = vec![SUB, infile, "-c", "1"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().success();
    }

    #[test]
    fn bam_no_index_fails_fetch() {
        let infile = "tests/cases/no_index.bam";
        let passed_args = vec![SUB, infile, "-c", "1", "--strategy", "fetch"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().failure();
    }

    #[test]
    fn bam_with_no_start_or_end_regions_and_missing_chromosomes() {
        let infile = "tests/cases/no_start_end.bam";
        let passed_args = vec![SUB, infile, "-c", "1"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().success();
    }

    #[test]
    fn test_make_program_id_unique_no_program() {
        let raw_header = "@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chromosome\tLN:5399960
@PG\tID:minimap2\tPN:minimap2\tVN:2.26-r1175\tCL:minimap2 -aL --cs --MD -t 4 -x map-ont KPC2__202310.5x.fq.gz
@PG\tID:samtools\tPN:samtools\tPP:minimap2\tVN:1.19.2\tCL:samtools sort -@ 4 -o KPC2__202310.5x.bam
@PG\tID:samtools.1\tPN:samtools\tPP:samtools\tVN:1.19\tCL:samtools view -s 0.5 -o test.bam KPC2__202310.5x.bam";
        let header = raw_header.parse().unwrap();
        let program_id = "rasusa";
        let actual = make_program_id_unique(&header, program_id);
        let expected = (
            Cow::<str>::Borrowed(program_id),
            Some("samtools.1".to_string()),
        );
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_make_program_id_unique_one_program_occurrence() {
        let raw_header = "@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chromosome\tLN:5399960
@PG\tID:minimap2\tPN:minimap2\tVN:2.26-r1175\tCL:minimap2 -aL --cs --MD -t 4 -x map-ont KPC2__202310.5x.fq.gz
@PG\tID:samtools\tPN:samtools\tPP:minimap2\tVN:1.19.2\tCL:samtools sort -@ 4 -o KPC2__202310.5x.bam
@PG\tID:samtools.1\tPN:samtools\tPP:samtools\tVN:1.19\tCL:samtools view -s 0.5 -o test.bam KPC2__202310.5x.bam";
        let header = raw_header.parse().unwrap();
        let program_id = "minimap2";
        let actual = make_program_id_unique(&header, program_id);
        let expected = (
            Cow::<str>::Owned("minimap2.1".to_string()),
            Some("samtools.1".to_string()),
        );
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_make_program_id_unique_two_program_occurrences() {
        let raw_header = "@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chromosome\tLN:5399960
@PG\tID:minimap2\tPN:minimap2\tVN:2.26-r1175\tCL:minimap2 -aL --cs --MD -t 4 -x map-ont KPC2__202310.5x.fq.gz
@PG\tID:samtools\tPN:samtools\tPP:minimap2\tVN:1.19.2\tCL:samtools sort -@ 4 -o KPC2__202310.5x.bam
@PG\tID:samtools.1\tPN:samtools\tPP:samtools\tVN:1.19\tCL:samtools view -s 0.5 -o test.bam KPC2__202310.5x.bam";

        let header = raw_header.parse().unwrap();
        let program_id = "samtools";
        let actual = make_program_id_unique(&header, program_id);
        let expected = (
            Cow::<str>::Owned("samtools.2".to_string()),
            Some("samtools.1".to_string()),
        );
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_make_program_id_unique_no_programs() {
        let raw_header = "@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chromosome\tLN:5399960";
        let header = raw_header.parse().unwrap();
        let program_id = "samtools";
        let actual = make_program_id_unique(&header, program_id);
        let expected = (Cow::Borrowed("samtools"), None);
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_original_issue_random_compare_insufficient_shuffling() {
        // This test demonstrates the original issue: random_compare doesn't properly shuffle
        // groups of records with the same position because it only introduces randomness
        // at the comparison level, not at the group level.
        // This relates to issue #76.

        use std::cmp::Ordering;

        // Simulate the old random_compare function
        fn old_random_compare<T: Ord>(a: T, b: T, rng: &mut impl Rng) -> Ordering {
            if a == b {
                // Introduce randomness when elements are equal
                if rng.gen::<bool>() {
                    Ordering::Less
                } else {
                    Ordering::Greater
                }
            } else {
                a.cmp(&b)
            }
        }

        // Simulate the old random_sort function
        fn old_random_sort<T, K: Ord + Copy>(
            vec: &mut [T],
            key_extractor: fn(&T) -> K,
            mut rng: impl Rng,
        ) {
            vec.sort_by(|a, b| old_random_compare(key_extractor(a), key_extractor(b), &mut rng));
        }

        let mut rng = StdRng::seed_from_u64(42);

        // Create records with the same position
        let mut record1 = RecordBuf::default();
        *record1.alignment_start_mut() = Position::new(100);
        *record1.name_mut() = Some("read1".parse().unwrap());

        let mut record2 = RecordBuf::default();
        *record2.alignment_start_mut() = Position::new(100);
        *record2.name_mut() = Some("read2".parse().unwrap());

        let mut record3 = RecordBuf::default();
        *record3.alignment_start_mut() = Position::new(100);
        *record3.name_mut() = Some("read3".parse().unwrap());

        let records = vec![&record1, &record2, &record3];
        let original_order: Vec<Vec<u8>> =
            records.iter().map(|r| r.name().unwrap().to_vec()).collect();

        // Test the old approach - it often fails to properly shuffle
        let mut same_order_count = 0;
        for _ in 0..20 {
            let mut test_records = records.clone();
            old_random_sort(
                &mut test_records,
                |record| record.alignment_start(),
                &mut rng,
            );

            let new_order: Vec<Vec<u8>> = test_records
                .iter()
                .map(|r| r.name().unwrap().to_vec())
                .collect();
            if new_order == original_order {
                same_order_count += 1;
            }
        }

        // The old approach often keeps the same order because random_compare
        // doesn't guarantee proper shuffling of equal elements
        println!("Old approach: {same_order_count} out of 20 iterations kept the same order");

        // Now test our new approach
        let mut new_same_order_count = 0;
        for _ in 0..20 {
            let mut test_records = records.clone();
            shuffle_records_by_position(&mut test_records, &mut rng);

            let new_order: Vec<Vec<u8>> = test_records
                .iter()
                .map(|r| r.name().unwrap().to_vec())
                .collect();
            if new_order == original_order {
                new_same_order_count += 1;
            }
        }

        println!("New approach: {new_same_order_count} out of 20 iterations kept the same order");

        // The new approach should shuffle much more effectively
        // (though due to randomness, it might occasionally keep the same order)
        assert!(
            new_same_order_count < same_order_count,
            "New shuffling approach should be more effective than old random_compare approach"
        );
    }

    #[test]
    fn test_make_program_id_unique_program_id_startswith_same_substring() {
        let raw_header = "@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chromosome\tLN:5399960
@PG\tID:minimap2\tPN:minimap2\tVN:2.26-r1175\tCL:minimap2 -aL --cs --MD -t 4 -x map-ont KPC2__202310.5x.fq.gz
@PG\tID:samtoolsfoo\tPN:samtools\tPP:minimap2\tVN:1.19.2\tCL:samtools sort -@ 4 -o KPC2__202310.5x.bam
@PG\tID:samtools\tPN:samtools\tPP:minimap2\tVN:1.19.2\tCL:samtools sort -@ 4 -o KPC2__202310.5x.bam
@PG\tID:samtoolsfoo.1\tPN:samtools\tPP:samtools\tVN:1.19\tCL:samtools view -s 0.5 -o test.bam KPC2__202310.5x.bam";
        let header = raw_header.parse().unwrap();
        let program_id = "samtools";
        let actual = make_program_id_unique(&header, program_id);
        let expected = (
            Cow::<str>::Owned("samtools.1".to_string()),
            Some("samtoolsfoo.1".to_string()),
        );
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_scored_read_ordering() {
        let r1 = RecordBuf::default();
        let small = ScoredRead {
            record: r1,
            key: 10,
            end: 1000,
        };

        let r2 = RecordBuf::default();
        let big = ScoredRead {
            record: r2,
            key: 99,
            end: 1000,
        };

        assert!(big > small, "Higher key should be 'Greater' in ordering");
        assert!(small < big, "Lower key should be 'Less' in ordering");
    }

    #[test]
    fn test_scored_read_ordering_tie() {
        let mut r1 = RecordBuf::default();
        *r1.name_mut() = Some("readA".parse().unwrap());
        let record1 = ScoredRead {
            record: r1,
            key: 21,
            end: 1000,
        };

        let mut r2 = RecordBuf::default();
        *r2.name_mut() = Some("readB".parse().unwrap());
        let record2 = ScoredRead {
            record: r2,
            key: 21,
            end: 1001,
        };

        assert!(record2 > record1);
    }

    #[test]
    fn test_output_never_exceeds_target_depth() {
        let input_path = PathBuf::from("tests/cases/test.bam");
        let target_depth = 3;

        let output_file = NamedTempFile::new().unwrap();
        let output_path = output_file.path().to_path_buf();

        let mut align = Alignment {
            aln: input_path,
            output: Some(output_path.clone()),
            output_type: Some(Format::Bam),
            coverage: target_depth,
            seed: Some(2109),
            strategy: SubsamplingStrategy::Stream,
            swap_distance: 5,
            step_size: 100, // it is not used
            batch_size: 10_000,
        };

        // run the sweep line algorithm
        align.run().expect("Subsampling failed");

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
    fn bam_is_not_sorted_fails() {
        let infile = "tests/cases/test_not_sorted.bam";
        let passed_args = vec![SUB, infile, "-c", "1"];
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME")).unwrap();

        cmd.args(passed_args).assert().failure();
    }

    // helper function to run subsamping aln and get the resulted read names
    fn run_aln_get_reads_result(
        input: &Path,
        seed: Option<u64>,
        strategy: SubsamplingStrategy,
    ) -> Vec<String> {
        let target_depth = 3;
        let out = NamedTempFile::new().unwrap();

        let mut aln1 = Alignment {
            aln: input.to_path_buf(),
            output: Some(out.path().to_path_buf()),
            output_type: Some(Format::Bam),
            coverage: target_depth,
            seed,
            strategy,
            swap_distance: 5,
            step_size: 100,
            batch_size: 10_000,
        };

        aln1.run().expect("Subsampling failed");

        let mut reader = alignment::io::reader::Builder::default()
            .build_from_path(out.path())
            .unwrap();
        let header = reader.read_header().unwrap();

        reader
            .records(&header)
            .map(|r| {
                let rec = r.unwrap();
                // Get name as String from the record
                String::from_utf8_lossy(rec.name().unwrap()).to_string()
            })
            .collect()
    }

    #[test]
    fn test_reproducibility_stream_same_seed() {
        let input_path = Path::new("tests/cases/test.bam");
        let seed = Some(2109);

        let names1 = run_aln_get_reads_result(input_path, seed, SubsamplingStrategy::Stream);
        let names2 = run_aln_get_reads_result(input_path, seed, SubsamplingStrategy::Stream);

        // comapre the length and the read names
        assert_eq!(names1.len(), names2.len(), "Different read count");
        assert_eq!(names1, names2, "Different reads selected")
    }

    #[test]
    fn test_reproducibility_fetch_same_seed() {
        let input_path = Path::new("tests/cases/test.bam");
        let seed = Some(2109);

        let names1 = run_aln_get_reads_result(input_path, seed, SubsamplingStrategy::Fetch);
        let names2 = run_aln_get_reads_result(input_path, seed, SubsamplingStrategy::Fetch);

        // comapre the length and the read names
        assert_eq!(names1.len(), names2.len(), "Different read count");
        assert_eq!(names1, names2, "Different reads selected")
    }

    #[test]
    fn test_reproducibility_stream_diff_seed() {
        let input_path = Path::new("tests/cases/test.bam");
        let seed1 = Some(21);
        let seed2 = Some(9);

        let names1 = run_aln_get_reads_result(input_path, seed1, SubsamplingStrategy::Stream);
        let names2 = run_aln_get_reads_result(input_path, seed2, SubsamplingStrategy::Stream);

        // comapre the length and the read names
        assert_ne!(names1, names2, "Same reads selected")
    }

    #[test]
    fn test_reproducibility_fetch_diff_seed() {
        let input_path = Path::new("tests/cases/test.bam");
        let seed1 = Some(21);
        let seed2 = Some(9);

        let names1 = run_aln_get_reads_result(input_path, seed1, SubsamplingStrategy::Fetch);
        let names2 = run_aln_get_reads_result(input_path, seed2, SubsamplingStrategy::Fetch);

        // comapre the length and the read names
        assert_ne!(names1, names2, "Same reads selected")
    }

    #[test]
    fn unknown_input_extension_fails() {
        let input = Path::new("tests/cases/test");

        let mut aln = Alignment {
            aln: input.to_path_buf(),
            output: None,
            output_type: Some(Format::Bam),
            coverage: 2,
            seed: Some(2109),
            strategy: SubsamplingStrategy::Stream,
            swap_distance: 5,
            step_size: 100, // not used
            batch_size: 10_000,
        };
        assert!(aln.run().is_err());
    }

    #[test]
    fn unknown_output_extension_fails() {
        let input = Path::new("tests/cases/test.bam");
        let output = Path::new("tests/cases/result");

        let mut aln = Alignment {
            aln: input.to_path_buf(),
            output: Some(output.to_path_buf()),
            output_type: None,
            coverage: 2,
            seed: Some(2109),
            strategy: SubsamplingStrategy::Stream,
            swap_distance: 5,
            step_size: 100, // not used
            batch_size: 10_000,
        };
        assert!(aln.run().is_err());
    }

    #[test]
    fn test_equality_record() {
        let mut r1 = RecordBuf::default();
        *r1.name_mut() = Some("readA".parse().unwrap());
        let record1 = ScoredRead {
            record: r1.clone(),
            key: 21,
            end: 1000,
        };

        let mut r2 = RecordBuf::default();
        *r2.name_mut() = Some("readA".parse().unwrap());
        let record2 = ScoredRead {
            record: r2.clone(),
            key: 21,
            end: 1000,
        };

        assert_eq!(record1, record2);
    }

    #[test]
    fn test_subsampling_statistics() {
        let input_path = PathBuf::from("tests/cases/test.bam");

        let target_depth = 3;

        let output_temp = NamedTempFile::new().unwrap();

        let mut align = Alignment {
            aln: input_path,
            output: Some(output_temp.path().to_path_buf()),
            output_type: Some(Format::Bam),
            coverage: target_depth,
            seed: Some(2109),
            strategy: SubsamplingStrategy::Stream,
            swap_distance: 5,
            step_size: 100, // not used
            batch_size: 10_000,
        };
        align.run().expect("Subsampling failed");

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

    #[test]
    fn test_paired_end_retention_stream() {
        let input_path = PathBuf::from("tests/cases/test.paired.bam");

        let target_depth: u32 = 2;
        let output = NamedTempFile::new().unwrap();

        // subsample to 2x depth
        let mut align = Alignment {
            aln: input_path,
            output: Some(output.path().to_path_buf()),
            output_type: Some(Format::Bam),
            coverage: target_depth,
            seed: Some(2109),
            strategy: SubsamplingStrategy::Stream,
            swap_distance: 5,
            step_size: 100,
            batch_size: 1000,
        };

        align.run().expect("Subsampling failed");

        // verify
        let mut reader = alignment::io::reader::Builder::default()
            .build_from_path(output.path())
            .unwrap();
        // read the header, get the length of chromosome
        let header = reader.read_header().unwrap();

        let mut r1_names: Vec<String> = Vec::new();
        let mut r2_names: Vec<String> = Vec::new();

        for result in reader.records(&header) {
            let record = result.unwrap();
            let name = record.name().unwrap().to_string();
            let flags = record.flags().unwrap();

            if flags.is_first_segment() {
                r1_names.push(name);
            } else if flags.is_last_segment() {
                r2_names.push(name);
            }
        }
        // make sure r1 and r2 names are identical because they come from the same template, sam spec guarantees this
        // read more: https://samtools.github.io/hts-specs/SAMv1.pdf
        r1_names.sort();
        r2_names.sort();
        assert_eq!(r1_names, r2_names, "Mismatch!");
    }

    #[test]
    fn test_paired_end_retention_fetch() {
        let input_path = PathBuf::from("tests/cases/test.paired.bam");

        let target_depth: u32 = 2;
        let output = NamedTempFile::new().unwrap();

        // subsample to 2x depth
        let mut align = Alignment {
            aln: input_path,
            output: Some(output.path().to_path_buf()),
            output_type: Some(Format::Bam),
            coverage: target_depth,
            seed: Some(2109),
            strategy: SubsamplingStrategy::Fetch,
            swap_distance: 5,
            step_size: 100,
            batch_size: 10000,
        };

        align.run().expect("Subsampling failed");

        // verify
        let mut reader = alignment::io::reader::Builder::default()
            .build_from_path(output.path())
            .unwrap();
        // read the header, get the length of chromosome
        let header = reader.read_header().unwrap();

        let mut r1_names: Vec<String> = Vec::new();
        let mut r2_names: Vec<String> = Vec::new();

        for result in reader.records(&header) {
            let record = result.unwrap();
            let name = record.name().unwrap().to_string();
            let flags = record.flags().unwrap();

            if flags.is_first_segment() {
                r1_names.push(name);
            } else if flags.is_last_segment() {
                r2_names.push(name);
            }
        }
        // make sure r1 and r2 names are identical because they come from the same template, sam spec guarantees this
        // read more: https://samtools.github.io/hts-specs/SAMv1.pdf
        r1_names.sort();
        r2_names.sort();

        assert_eq!(r1_names, r2_names, "Mismatch!");
    }
}
