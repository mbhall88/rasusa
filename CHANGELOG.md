# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- Switch from using `snafu` and `failure` for error handling to `anyhow` and `thiserror`. Based on the procedure outlined in [this excellent blog post][error-blog].

## [0.3.0]

Version 0.3.0 may give different results to previous versions. If so,
the differences will likely be a handful of extra reads (possibly none).
The reason for this is `--coverage` is now treated as a float.
Previously we immediately round coverage down to the nearest integer. As
the number of reads to keep is based on the target total number of
bases, which is coverage * genome size. So if coverage is 10.7 and
genome size is 100, previously our target number of bases would have
been 1000, whereas now, it would be 1070.

### Changed
- `--coverage` is now treated as a `f32` instead of being converted
  immediately to an integer [#19][19].
- Updated `rust-bio` to version 0.31.0. This means `rasusa` now handles
  wrapped fastq files.
- Preallocate fastx records instead of using iterator. Gives marginal
  speedup.
- Added `bash` to the docker image b47a8b75943098bdd845b7758cf2eab01ef5a3d8

## [0.2.0]

### Added
- Support paired Illumina [#15](https://github.com/mbhall88/rasusa/issues/15)


[unreleased]: https://github.com/mbhall88/rasusa/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/mbhall88/rasusa/releases/tag/v0.3.0
[0.2.0]: https://github.com/mbhall88/rasusa/releases/tag/v0.2.0
[19]: https://github.com/mbhall88/rasusa/issues/19
[error-blog]: https://nick.groenen.me/posts/rust-error-handling/