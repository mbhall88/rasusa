# Requirements: Support Unaligned SAM/BAM/CRAM in `reads` subcommand

## Overview
The `reads` subcommand currently only supports FASTA/FASTQ files via `needletail`. This requirement is to add support for unaligned SAM, BAM, and CRAM files as both input and output, using the `noodles` library.

## User Stories
- As a bioinformatician, I want to subsample my unaligned BAM files (e.g., from Nanopore or PacBio) directly without converting to FASTQ first.
- As a user with paired-end unaligned BAM data, I want to subsample it while ensuring that both reads of a pair are either kept or discarded together.
- I want the subsampling to correctly account for the total number of bases or target coverage when using unaligned SAM/BAM/CRAM inputs.

## Functional Requirements
- **Format Detection**: Automatically detect if an input file is SAM, BAM, or CRAM.
- **Single-Pass/Two-Pass Support**: Maintain the current two-pass approach for the `reads` subcommand (first pass for lengths, second for filtering).
- **Paired-End Integrity**: If two unaligned SAM/BAM/CRAM files are provided as input (paired-end), ensure matching records are handled together.
- **Single-File Paired-End**: Support paired-end data stored in a single unaligned SAM/BAM/CRAM file (where pairs are consecutive or identified by name). *Note: Use the `SEGMENTED` flag (bit 0x1) for detection.*
- **Output Support**: Allow writing the subsampled output back to unaligned SAM/BAM/CRAM (matching input format or as specified).
- **Target Calculation**: Correctly calculate target bases/reads for these new formats.

## Non-Functional Requirements
- **Performance**: Minimize overhead when using `noodles` for unaligned records.
- **Memory**: Be mindful of memory usage, especially if storing read names for pairing.
- **Compatibility**: Ensure no regressions for existing FASTA/FASTQ support.

## Technical Constraints
- Use the `noodles` library for all unaligned SAM/BAM/CRAM I/O.
- Maintain consistency with the `Runner` trait and existing `reads` subcommand structure.
