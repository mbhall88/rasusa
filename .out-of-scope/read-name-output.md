# Writing the names of kept reads to a separate output

Rasusa does not have an option that writes the names (IDs) of the reads it keeps to a separate file or to the log.

## Why this is out of scope

Rasusa already provides two ways to recover the kept reads' names, so a dedicated output is redundant:

1. **The seed records the selection.** Given the same input, seed, and Rasusa version, the output is byte-identical. Re-running with `-s <seed>` reproduces exactly which reads were kept, so the seed is the record.
2. **The output stream already carries the names.** When streaming into another tool, `tee` the output and extract the names:

   ```sh
   hifiasm -o asm <(
     rasusa reads -p 0.5 -s 12345 -O fastq reads.bam \
       | tee >(awk 'NR % 4 == 1 { print substr($1, 2) }' > kept.reads.list)
   )
   ```

   (Rasusa before 5.1.1 also needs `-o /dev/stdout` here: see #238.)

A dedicated `--read-names <FILE>` flag would add CLI surface and would have to be threaded through every place Rasusa writes output (FASTX two-pass, one-pass single-end and paired, and alignment-source two-pass and one-pass). It would also need defined semantics for paired templates. For a need these workarounds already cover, that cost is not justified.

Emitting names through trace-level logging was also considered and rejected. Logs go to unbuffered stderr mixed with other messages, so it is slow at millions of reads and turns the log format into an accidental interface.

## Reconsidering

This is a demand question, not a design objection. If enough users need the names as a first-class output, the preferred design is a dedicated flag writing one name per **template** (not logging), consistent with ADR-0001's constant-memory constraint.

## Prior requests

- #237: "Add an option to write the names of subsampled reads to a separate file"
