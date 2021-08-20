# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.5.0]

### Added

- Support for LZMA, Bzip, and Gzip output compression (thanks to
  [`niffler`](https://github.com/luizirber/niffler/)). This is either inferred from the
  file extension or manually via the `-O` option.
- Option to specify the compression level for the output via `-l`

### Changed

- Use a `Vec<bool>` instead of `HashSet` to store the indices of reads to keep. This
  gives a nice little speedup (see [#28][28]), A big thank you to
  [@natir](https://github.com/natir) for this.

### Fixed
- Restore compression of output files [[#27][27]]

## [0.4.2]

### Fixed

- I had stupidly forgetten to merge the fix for [#22][22] onto master 🤦

## [0.4.1]

### Fixes

- Releasing cross-compiled binaries didn't work for version 0.4.0
- Docker image is now correctly built

## [0.4.0]

### Changed

- Switch from using `snafu` and `failure` for error handling to `anyhow` and
  `thiserror`. Based on the procedure outlined in [this excellent blog
  post][error-blog].
- Switched fasta/q parsing to use [needletail](https://github.com/onecodex/needletail)
  instead of rust-bio. See [benchmark] for improvement in runtimes.
- Changed the way Illumina paired reads are subsampled. Previously, there was an
  assumption made that the reads of a pair were both the same length as the R1 read. We
  are now more careful and look at each read's length individually [[#22][22]]
- Moved container hosting to quay.io

## [0.3.0]

Version 0.3.0 may give different results to previous versions. If so, the differences
will likely be a handful of extra reads (possibly none). The reason for this is
`--coverage` is now treated as a float. Previously we immediately round coverage down to
the nearest integer. As the number of reads to keep is based on the target total number
of bases, which is coverage * genome size. So if coverage is 10.7 and genome size is
100, previously our target number of bases would have been 1000, whereas now, it would
be 1070.

### Changed

- `--coverage` is now treated as a `f32` instead of being converted immediately to an
  integer [#19][19].
- Updated `rust-bio` to version 0.31.0. This means `rasusa` now handles wrapped fastq
  files.
- Preallocate fastx records instead of using iterator. Gives marginal speedup.
- Added `bash` to the docker image b47a8b75943098bdd845b7758cf2eab01ef5a3d8

## [0.2.0]

### Added

- Support paired Illumina [#15](https://github.com/mbhall88/rasusa/issues/15)

[0.2.0]: https://github.com/mbhall88/rasusa/releases/tag/0.2.0
[0.3.0]: https://github.com/mbhall88/rasusa/releases/tag/0.3.0
[0.4.0]: https://github.com/mbhall88/rasusa/releases/tag/0.4.0
[0.4.1]: https://github.com/mbhall88/rasusa/releases/tag/0.4.1
[0.4.2]: https://github.com/mbhall88/rasusa/releases/tag/0.4.2
[0.5.0]: https://github.com/mbhall88/rasusa/releases/tag/0.5.0
[19]: https://github.com/mbhall88/rasusa/issues/19
[22]: https://github.com/mbhall88/rasusa/issues/22
[27]: https://github.com/mbhall88/rasusa/issues/27
[28]: https://github.com/mbhall88/rasusa/pull/28
[benchmark]: https://github.com/mbhall88/rasusa#benchmark
[error-blog]: https://nick.groenen.me/posts/rust-error-handling/
[unreleased]: https://github.com/mbhall88/rasusa/compare/0.5.0...HEAD

