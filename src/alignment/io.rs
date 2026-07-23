use std::fs::File;
use std::io::{self, Write};

use anyhow::{Context, Result};
use log::info;
use rand::prelude::*;
use rand::random;

use noodles::sam::Header;

use crate::format::infer_format_from_path;
use crate::threading::build_alignment_writer;

use super::args::Alignment;
use super::header::program_entry;

// re-exported so `stream.rs`/`fetch.rs`/`mod.rs` can keep referring to it as `super::io::AlignmentWriter`
pub(super) use crate::threading::AlignmentWriter;

impl Alignment {
    pub(super) fn setup_resources(
        &self,
        input_header: &Header,
    ) -> Result<(rand_pcg::Pcg64, AlignmentWriter)> {
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
        let (pg_id, pg_map) = program_entry(&header);

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

        let output_fmt = match &self.output_format {
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

        // use Box<dyn Write + Send> to make File and Stdout compatible, and to allow wrapping in
        // a multithreaded BGZF encoder when `--threads` > 1.
        let sink: Box<dyn Write + Send> = match &self.output {
            Some(path) => {
                let path = path.as_path();
                info!("Writing subsampled alignment to: {:?}", path);
                let file = File::create(path).context("Failed to create output alignment file")?;
                Box::new(file)
            }
            None => {
                info!("Writing subsampled alignment to stdout");
                // `Stdout` (rather than a `StdoutLock`) so the sink is `Send`: it can be wrapped
                // in a multithreaded BGZF encoder, whose background writer thread needs to own
                // it. `Stdout::write` still locks internally on every call.
                Box::new(io::stdout())
            }
        };

        let mut writer = build_alignment_writer(sink, output_fmt, self.threads)?;

        writer.write_header(&header)?;

        Ok((rng, writer))
    }
}
