![rasusa](img/logo.png)

[![Build Status](https://github.com/mbhall88/rasusa/actions/workflows/rust-ci.yaml/badge.svg?branch=main)](https://github.com/mbhall88/rasusa/actions/workflows/rust-ci.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![github release version](https://img.shields.io/github/v/release/mbhall88/rasusa)](https://github.com/mbhall88/rasusa/releases)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03941/status.svg)](https://doi.org/10.21105/joss.03941)

**Ra**ndomly **su**b**sa**mple sequencing reads or alignments.

> Hall, M. B., (2022). Rasusa: Randomly subsample sequencing reads to a specified coverage. Journal of Open Source
> Software, 7(69), 3941, https://doi.org/10.21105/joss.03941

[TOC]: #

## Table of Contents

- [Table of Contents](#table-of-contents)
- [Motivation](#motivation)
- [Install](#install)
  - [`cargo`](#cargo)
  - [`conda`](#conda)
  - [Container](#container)
    - [`singularity`](#singularity)
    - [`docker`](#docker)
  - [`homebrew`](#homebrew)
  - [Build locally](#build-locally)
- [Usage](#usage)
  - [Basic usage - reads](#basic-usage---reads)
  - [Basic usage - alignments](#basic-usage---alignments)
  - [Required parameters](#required-parameters)
    - [Input](#input)
    - [Coverage](#coverage)
      - [`-c`, `--coverage`](#-c---coverage)
    - [Genome size](#genome-size)
      - [`-g`, `--genome-size`](#-g---genome-size)
  - [Optional parameters](#optional-parameters)
    - [Output](#output)
      - [`-o`, `--output`](#-o---output)
    - [Output compression/format](#output-compressionformat)
      - [`-O`, `--output-type`](#-o---output-type)
    - [Compresion level](#compresion-level)
      - [`-l`, `--compress-level`](#-l---compress-level)
    - [Target number of bases](#target-number-of-bases)
      - [`-b`, `--bases`](#-b---bases)
    - [Number of reads](#number-of-reads)
      - [`-n`, `--num`](#-n---num)
    - [Fraction of reads](#fraction-of-reads)
      - [`-f`, `--frac`](#-f---frac)
    - [Random seed](#random-seed)
      - [`-s`, `--seed`](#-s---seed)
    - [Verbosity](#verbosity)
      - [`-v`](#-v)
  - [Full usage](#full-usage)
    - [`reads` command](#reads-command)
    - [`aln` command](#aln-command)
- [Benchmark](#benchmark)
  - [Single long read input](#single-long-read-input)
    - [Results](#results)
  - [Paired-end input](#paired-end-input)
    - [Results](#results-1)
- [Contributing](#contributing)
- [Citing](#citing)
  - [Bibtex](#bibtex)

## Motivation

I couldn't find a tool for subsampling fastq reads that met my requirements. All the
strategies I could find fell short as they either just wanted a number or percentage of
reads to subsample to or, if they did subsample to a coverage, they assume all reads are
the same size (i.e Illumina). As I mostly work with long-read data this posed a problem
if I wanted to subsample a file to certain coverage, as length of reads was never taken
into account. `rasusa` addresses this shortcoming.

A workaround I had been using for a while was using [`filtlong`][filtlong]. It was
simple enough, I just figure out the number of bases I need to achieve a (theoretical)
coverage for my sample. Say I have a fastq from an _E. coli_ sample with 5 million reads
and I want to subset it to 50x coverage. I just need to multiply the expected size of
the sample's genome, 4.6 million base pairs, by the coverage I want and I have my target
bases - 230 million base pairs. In `filtlong`, I can do the following

```sh
target=230000000
filtlong --target_bases "$target" reads.fq > reads.50x.fq
```

However, this is technically not the intended function of `filtlong`; it's a quality
filtering tool. What you get in the end is a subset of the ["highest scoring"][score]
reads at a (theoretical) coverage of 50x. Depending on your circumstances, this might be
what you want. However, you bias yourself towards the best/longest reads in the dataset
\- not a fair representation of your dataset as a whole. There is also the possibility
of favouring regions of the genome that produce longer/higher quality reads. [De Maio
_et al._][mgen-ref] even found that by randomly subsampling nanopore reads you achieve
_better_ genome assemblies than if you had filtered.

So, depending on your circumstances, an unbiased subsample of your reads might be what
you need. And if this is the case, `rasusa` has you covered.

## Install

![GitHub Downloads (all assets, all releases)](https://img.shields.io/github/downloads/mbhall88/rasusa/total)

**tl;dr**: precompiled binary

```shell
curl -sSL rasusa.mbh.sh | sh
# or with wget
wget -nv -O - rasusa.mbh.sh | sh
```

You can also pass options to the script like so

```
$ curl -sSL rasusa.mbh.sh | sh -s -- --help
install.sh [option]

Fetch and install the latest version of rasusa, if rasusa is already
installed it will be updated to the latest version.

Options
        -V, --verbose
                Enable verbose output for the installer

        -f, -y, --force, --yes
                Skip the confirmation prompt during installation

        -p, --platform
                Override the platform identified by the installer [default: apple-darwin]

        -b, --bin-dir
                Override the bin installation directory [default: /usr/local/bin]

        -a, --arch
                Override the architecture identified by the installer [default: x86_64]

        -B, --base-url
                Override the base URL used for downloading releases [default: https://github.com/mbhall88/rasusa/releases]

        -h, --help
                Display this help message
```

### `cargo`

[![Crates.io](https://img.shields.io/crates/v/rasusa.svg)](https://crates.io/crates/rasusa)
![Crates.io Total Downloads](https://img.shields.io/crates/d/rasusa)

Prerequisite: [`rust` toolchain][rust] (min. v1.74.1)

```sh
cargo install rasusa
```

### `conda`

[![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/rasusa)](https://anaconda.org/bioconda/rasusa)
[![bioconda version](https://anaconda.org/bioconda/rasusa/badges/platforms.svg)](https://anaconda.org/bioconda/rasusa)
![Conda](https://img.shields.io/conda/dn/bioconda/rasusa)

Prerequisite: [`conda`][conda] (and bioconda channel [correctly set up][channels])

```sh
conda install rasusa
```

Thank you to Devon Ryan ([@dpryan79][dpryan79]) for [help debugging the bioconda
recipe][pr-help].

### Container

Docker images are hosted at [quay.io]. For versions 0.3.0 and earlier, the images were
hosted on [Dockerhub][dockerhub].

#### `singularity`

Prerequisite: [`singularity`][singularity]

```sh
URI="docker://quay.io/mbhall88/rasusa"
singularity exec "$URI" rasusa --help
```

The above will use the latest version. If you want to specify a version then use a
[tag][quay.io] (or commit) like so.

```sh
VERSION="0.8.0"
URI="docker://quay.io/mbhall88/rasusa:${VERSION}"
```

#### `docker`

[![Docker Repository on Quay](https://quay.io/repository/mbhall88/rasusa/status "Docker Repository on Quay")](https://quay.io/repository/mbhall88/rasusa)

Prerequisite: [`docker`][docker]

```sh
docker pull quay.io/mbhall88/rasusa
docker run quay.io/mbhall88/rasusa --help
```

You can find all the available tags on the [quay.io repository][quay.io]. Note: versions
prior to 0.4.0 were housed on [Docker Hub](https://hub.docker.com/r/mbhall88/rasusa).

### `homebrew`

![Homebrew](https://img.shields.io/homebrew/v/rasusa)
![Homebrew Downloads](https://img.shields.io/homebrew/installs/dy/rasusa)

Prerequisite: [`homebrew`][homebrew]

```sh
brew install rasusa
```

### Build locally

Prerequisite: [`rust` toolchain][rust]

```sh
git clone https://github.com/mbhall88/rasusa.git
cd rasusa
cargo build --release
target/release/rasusa --help
# if you want to check everything is working ok
cargo test --all
```

## Usage

### Basic usage - reads

Subsample fastq reads

```
rasusa reads --coverage 30 --genome-size 4.6mb in.fq
```

The above command will output the subsampled file to `stdout`.

Or, if you have paired Illumina

```
rasusa reads --coverage 30 --genome-size 4g -o out.r1.fq -o out.r2.fq r1.fq r2.fq
```

For more details on the above options, and additional options, see below.

### Basic usage - alignments

Subsample alignments

```
rasusa aln --coverage 30 in.bam | samtools sort -o out.bam
```

this will subsample each position in the alignment to 30x coverage.

### Required parameters

There are three required options to run `rasusa reads`.

#### Input

This positional argument specifies the file(s) containing the reads or alignments you would like to subsample. The
file(s) must be valid fasta or fastq format for the `reads` command and can be compressed (with a tool such as
`gzip`). For the `aln` command, the file must be a valid **indexed** SAM/BAM file.  
If two files are passed to `reads`, `rasusa` will assume they are paired-end reads.

> Bash wizard tip 🧙: Let globs do the work for you `r*.fq`

#### Coverage

##### `-c`, `--coverage`

> Not required if [`--bases`](#target-number-of-bases) is present for `reads`

This option is used to determine the minimum coverage to subsample the reads to. For the `reads` command, it can
be specified as an integer (100), a decimal/float (100.7), or either of the previous
suffixed with an 'x' (100x). For the `aln` command, it is an integer only.

Due to the method for determining how many bases are required to achieve the
desired coverage in the `reads` command, the actual coverage, in the end, could be slightly higher than
requested. For example, if the last included read is very long. The log messages should
inform you of the actual coverage in the end.

For the `aln` command, the coverage is the minimum number of reads that should be present at each position in the
alignment. If a position has fewer than the requested number of reads, all reads at that position will be included. In
addition, there will be (small) regions with more than the requested number of reads - usually localised to where the
alignment of a read ends. This is
because when the alignment of a selected read ends, the next read is selected based on it spanning the end of the
previous alignment.
When selecting this next alignment, we preference alignments whose start is closest to the end of the previous
alignment, ensuring minimal overlap with the previous alignment. See the below screenshot from IGV for a visual example.

![IGV screenshot](img/igv_panel.png)

#### Genome size

##### `-g`, `--genome-size`

> Not valid for `aln`

> Not required if [`--bases`](#target-number-of-bases) is present for `reads`

The genome size of the input is also required. It is used to determine how many bases
are necessary to achieve the desired coverage. This can, of course, be as precise or
rough as you like.  
Genome size can be passed in many ways. As a plain old integer (1600), or with a metric
suffix (1.6kb). All metric suffixes can have an optional 'b' suffix and be lower, upper,
or mixed case. So 'Kb', 'kb' and 'k' would all be inferred as 'kilo'. Valid metric
suffixes include:

- Base (b) - multiplies by 1
- Kilo (k) - multiplies by 1,000
- Mega (m) - multiplies by 1,000,000
- Giga (g) - multiplies by 1,000,000,000
- Tera (t) - multiplies by 1,000,000,000,000

Alternatively, a [FASTA/Q index file][faidx] can be given and the genome size will be
set to the sum of all reference sequences in it.

> [!TIP]
> If you want to use `rasusa` in a scenario where you don't know what the genome size is, 
> such as in an automated pipeline that can take in any kind of organism, you could estimate 
> the genome size with something like [`lrge`](https://github.com/mbhall88/lrge) (#shamelessplug).
> 
> ```
> $ gsize=$(lrge reads.fq)
> $ rasusa reads -g $gsize -c 10 reads.fq
> ```
> `lrge` is designed for long reads. If you want to estimate the genome size from short 
> reads, you could use something like [Mash](https://github.com/marbl/Mash) or 
> [GenomeScope2](https://github.com/tbenavi1/genomescope2.0). See [the `lrge` docs](https://github.com/tbenavi1/genomescope2.0) 
> for examples of how Mash/GenomeScope2 can be used for this task.

[faidx]: https://www.htslib.org/doc/faidx.html

### Optional parameters

#### Output

##### `-o`, `--output`

**`reads`**

NOTE: This parameter is required if passing paired Illumina data to `reads`.

By default, `rasusa` will output the subsampled file to `stdout` (if one file is given).
If you would prefer to specify an output file path, then use this option.

Output for Illumina paired files must be specified using `--output` twice - `-o out.r1.fq -o out.r2.fq`

The ordering of the output files is assumed to be the same as the input.  
_Note: The output will always be in the same format as the input. You cannot pass fastq
as input and ask for fasta as output._

`rasusa reads` will also attempt to automatically infer whether compression of the output
file(s) is required. It does this by detecting any of the supported extensions:

- `.gz`: will compress the output with [`gzip`][gzip]
- `.bz` or `.bz2`: will compress the output with [`bzip2`][bzip]
- `.lzma`: will compress the output with the [`xz`][xz] LZMA algorithm

**`aln`**

For the `aln` command, the output file format will be the same as the input if writing to stdout, otherwise it will be
inferred from the file extension. 

**Note:** the output alignment will most likely **not be sorted**. You can use `samtools sort` to sort the output. e.g.,

```
rasusa aln -c 5 in.bam | samtools sort -o out.bam
```

[gzip]: http://www.gzip.org/

[bzip]: https://sourceware.org/bzip2/

[xz]: https://tukaani.org/xz/

#### Output compression/format

##### `-O`, `--output-type`

**`reads`**

Use this option to manually set the compression algoritm to use for the output file(s).
It will override any format automatically detected from the output path.

Valid options are:

- `g`: [`gzip`][gzip]
- `b`: [`bzip2`][bzip]
- `l` or `x`: [`xz`][xz] LZMA algorithm
- `z`: [`zstd`][zstd]
- `u`: no compression

**`aln`**

Use this option to manually set the output file format. By default, the same format as the input will be used, or the
format will be guessed from the `--output` path extension if given. Valid options are:

- `b` or `bam`: BAM
- `c` or `cram`: CRAM
- `s` or `sam`: SAM

*Note: all values to this option are case insensitive.*

#### Compresion level

##### `-l`, `--compress-level`

Compression level to use if compressing the output. By default this is set to the default for the compression type being
output.

#### Target number of bases

##### `-b`, `--bases`

> `reads` only

Explicitly set the number of bases required in the subsample. This option takes the
number in the same format as [genome size](#genome-size).

*Note: if this option is given, genome size and coverage are not required, or ignored if
they are provided.*

#### Number of reads

##### `-n`, `--num`

> `reads` only

Explicitly set the number of reads in the subsample. This option takes the number in
the same format as [genome size](#genome-size).

When providing paired reads as input, this option will sample this many total read
pairs. For example, when passing `-n 20 r1.fq r2.fq`, the two output files will have
20 reads each, and the read ids will be the same in both.

*Note: if this option is given, genome size and coverage are not required.*

#### Fraction of reads

##### `-f`, `--frac`

> `reads` only

Explicitly set the fraction of total reads in the subsample. The value given to this
option can be a float or a percentage - i.e., `-f 0.5` and `-f 50` will both take half
of the reads.

*Note: if this option is given, genome size and coverage are not required.*

#### Random seed

##### `-s`, `--seed`

This option allows you to specify the [random seed][seed] used by the random subsampler. By explicitly setting this 
parameter, you make the subsample for the input reproducible. You only need to pass this parameter if you are likely 
to want to subsample the same input file again in the future and want the same subset of reads. However, if you forget 
to use this option, the seed generated by the system will be printed to the log output, allowing you to use it in the
future.

#### Strictly adhere to requested coverage

##### `-e`, `--strict`

> `reads` only

If the requested coverage, total bases, number of reads, or fraction of reads cannot be met, an error will be thrown. 
By default, a warning is displayed, and the maximum possible coverage, total bases, number of reads, or fraction of reads is used.

#### Verbosity

##### `-v`

Adding this optional flag will make the logging more verbose. By default, logging will
produce messages considered "info" or above (see [here][log-lvl] for more details). If
verbosity is switched on, you will additionally get "debug" level logging messages.

### Full usage

```text
$ rasusa --help
Randomly subsample reads or alignments

Usage: rasusa [OPTIONS] <COMMAND>

Commands:
  reads  Randomly subsample reads
  aln    Randomly subsample alignments to a specified depth of coverage
  cite   Get a bibtex formatted citation for this package
  help   Print this message or the help of the given subcommand(s)

Options:
  -v             Switch on verbosity
  -h, --help     Print help
  -V, --version  Print version
```

#### `reads` command

```text
$ rasusa reads --help
Randomly subsample reads

Usage: rasusa reads [OPTIONS] <FILE(S)>...

Arguments:
  <FILE(S)>...
          The fast{a,q} file(s) to subsample.

          For paired Illumina, the order matters. i.e., R1 then R2.

Options:
  -o, --output <OUTPUT>
          Output filepath(s); stdout if not present.

          For paired Illumina pass this flag twice `-o o1.fq -o o2.fq`

          NOTE: The order of the pairs is assumed to be the same as the input - e.g., R1 then R2. This option is required for paired input.
          
  -g, --genome-size <size|faidx>
          Genome size to calculate coverage with respect to. e.g., 4.3kb, 7Tb, 9000, 4.1MB

          Alternatively, a FASTA/Q index file can be provided and the genome size will be set to the sum of all reference sequences.

          If --bases is not provided, this option and --coverage are required

  -c, --coverage <FLOAT>
          The desired depth of coverage to subsample the reads to

          If --bases is not provided, this option and --genome-size are required

  -b, --bases <bases>
          Explicitly set the number of bases required e.g., 4.3kb, 7Tb, 9000, 4.1MB

          If this option is given, --coverage and --genome-size are ignored

  -n, --num <INT>
          Subsample to a specific number of reads

          If paired-end reads are passed, this is the number of (matched) reads from EACH file. This option accepts the same format as genome size - e.g., 1k will take 1000 reads

  -f, --frac <FLOAT>
          Subsample to a fraction of the reads - e.g., 0.5 samples half the reads

          Values >1 and <=100 will be automatically converted - e.g., 25 => 0.25

  -e, --strict
          Exit with an error if the requested coverage/bases/reads is not possible

  -s, --seed <INT>
          Random seed to use

  -v
          Switch on verbosity

  -O, --output-type <u|b|g|l|x|z>
          u: uncompressed; b: Bzip2; g: Gzip; l: Lzma; x: Xz (Lzma); z: Zstd

          Rasusa will attempt to infer the output compression format automatically from the filename extension. This option is used to override that. If writing to stdout, the default is uncompressed

  -l, --compress-level <1-21>
          Compression level to use if compressing output. Uses the default level for the format if not specified

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
```

#### `aln` command

```text
$ rasusa aln --help
Randomly subsample alignments to a specified depth of coverage

Usage: rasusa aln [OPTIONS] --coverage <INT> <FILE>

Arguments:
  <FILE>
          Path to the indexed alignment file (SAM/BAM/CRAM) to subsample

Options:
  -o, --output <FILE>
          Path to the output subsampled alignment file. Defaults to stdout (same format as input)

          The output is not guaranteed to be sorted. We recommend piping the output to `samtools sort`

  -O, --output-type <FMT>
          Output format. Rasusa will attempt to infer the format from the output file extension if not provided

  -c, --coverage <INT>
          The desired depth of coverage to subsample the alignment to

  -s, --seed <INT>
          Random seed to use

      --step-size <INT>
          When a region has less than the desired coverage, the step size to move along the chromosome to find more reads.

          The lowest of the step and the minimum end coordinate of the reads in the region will be used. This parameter can have a significant impact on the runtime of the subsampling process.

          [default: 100]

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
```

## Benchmark

> “Time flies like an arrow; fruit flies like a banana.”  
> ― Anthony G. Oettinger

The real question is: will `rasusa` just needlessly eat away at your precious time on
earth?

To do this benchmark, I am going to use [hyperfine][hyperfine].

The data I used comes from

> [Bainomugisa, Arnold, et al. "A complete high-quality MinION nanopore assembly of an
> extensively drug-resistant Mycobacterium tuberculosis Beijing lineage strain
> identifies novel variation in repetitive PE/PPE gene regions." Microbial genomics 4.7
> (2018).][1]

_Note, these benchmarks are for `reads` only as there is no other tool that replicates the functionality of `aln`._

### Single long read input

Download and rename the fastq

```shell
URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR649/008/SRR6490088/SRR6490088_1.fastq.gz"
wget "$URL" -O - | gzip -d -c > tb.fq
```

The file size is 2.9G, and it has 379,547 reads.  
We benchmark against `filtlong` using the same strategy outlined in
[Motivation](#motivation).

```shell
TB_GENOME_SIZE=4411532
COVG=50
TARGET_BASES=$(( TB_GENOME_SIZE * COVG ))
FILTLONG_CMD="filtlong --target_bases $TARGET_BASES tb.fq"
RASUSA_CMD="rasusa reads tb.fq -c $COVG -g $TB_GENOME_SIZE -s 1"
hyperfine --warmup 3 --runs 10 --export-markdown results-single.md \
     "$FILTLONG_CMD" "$RASUSA_CMD" 
```

#### Results

| Command                                      |       Mean [s] | Min [s] | Max [s] |     Relative |
|:---------------------------------------------|---------------:|--------:|--------:|-------------:|
| `filtlong --target_bases 220576600 tb.fq`    | 21.685 ± 0.055 |  21.622 |  21.787 | 21.77 ± 0.29 |
| `rasusa reads tb.fq -c 50 -g 4411532 -s 1` | 0.996 ±  0.013 |   0.983 |   1.023 |         1.00 |

**Summary**: `rasusa` ran 21.77 ± 0.29 times faster than `filtlong`.

### Paired-end input

Download and then deinterleave the fastq with [`pyfastaq`][pyfastaq]

```shell
URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/008/SRR6488968/SRR6488968.fastq.gz"
wget "$URL" -O - | gzip -d -c - | fastaq deinterleave - r1.fq r2.fq
```

Each file's size is 179M and has 283,590 reads.  
For this benchmark, we will use [`seqtk`][seqtk]. We will also test `seqtk`'s 2-pass
mode as this is analogous to `rasusa reads`.

```shell
NUM_READS=140000
SEQTK_CMD_1="seqtk sample -s 1 r1.fq $NUM_READS > /tmp/r1.fq; seqtk sample -s 1 r2.fq $NUM_READS > /tmp/r2.fq;"
SEQTK_CMD_2="seqtk sample -2 -s 1 r1.fq $NUM_READS > /tmp/r1.fq; seqtk sample -2 -s 1 r2.fq $NUM_READS > /tmp/r2.fq;"
RASUSA_CMD="rasusa reads r1.fq r2.fq -n $NUM_READS -s 1 -o /tmp/r1.fq -o /tmp/r2.fq"
hyperfine --warmup 10 --runs 100 --export-markdown results-paired.md \
     "$SEQTK_CMD_1" "$SEQTK_CMD_2" "$RASUSA_CMD"
```

#### Results

| Command                                                                                           |     Mean [ms] | Min [ms] | Max [ms] |    Relative |
|:--------------------------------------------------------------------------------------------------|--------------:|---------:|---------:|------------:|
| `seqtk sample -s 1 r1.fq 140000 > /tmp/r1.fq; seqtk sample -s 1 r2.fq 140000 > /tmp/r2.fq;`       |  907.7 ± 23.6 |    875.4 |    997.8 | 1.84 ± 0.62 |
| `seqtk sample -2 -s 1 r1.fq 140000 > /tmp/r1.fq; seqtk sample -2 -s 1 r2.fq 140000 > /tmp/r2.fq;` |  870.8 ± 54.9 |    818.2 |   1219.8 | 1.77 ± 0.61 |
| `rasusa reads r1.fq r2.fq -n 140000 -s 1 -o /tmp/r1.fq -o /tmp/r2.fq`                           | 492.2 ± 165.4 |    327.4 |    887.4 |        1.00 |

**Summary**: `rasusa reads` ran 1.84 times faster than `seqtk` (1-pass) and 1.77 times faster
than `seqtk` (2-pass)

So, `rasusa reads` is faster than `seqtk` but doesn't require a fixed number of reads -
allowing you to avoid doing maths to determine how many reads you need to downsample to
a specific coverage. 🤓

## Contributing

If you would like to help improve `rasusa` you are very welcome!

For changes to be accepted, they must pass the CI and coverage checks. These include:

- Code is formatted with `rustfmt`. This can be done by running `cargo fmt` in the
  project directory.
- There are no compiler errors/warnings. You can check this by running `cargo clippy
  --all-features --all-targets -- -D warnings`
- Code coverage has not reduced. If you want to check coverage before pushing changes, I
  use [`kcov`][kcov].

## Citing

If you use `rasusa` in your research, it would be very much appreciated if you could
cite it.

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03941/status.svg)](https://doi.org/10.21105/joss.03941)

> Hall, M. B., (2022). Rasusa: Randomly subsample sequencing reads to a specified coverage. Journal of Open Source
> Software, 7(69), 3941, https://doi.org/10.21105/joss.03941

### Bibtex

You can get the following citation by running `rasusa cite`

```Bibtex
@article{Hall2022,
  doi = {10.21105/joss.03941},
  url = {https://doi.org/10.21105/joss.03941},
  year = {2022},
  publisher = {The Open Journal},
  volume = {7},
  number = {69},
  pages = {3941},
  author = {Michael B. Hall},
  title = {Rasusa: Randomly subsample sequencing reads to a specified coverage},
  journal = {Journal of Open Source Software}
}

```

[1]: https://doi.org/10.1099/mgen.0.000188

[brew-tap]: https://github.com/brewsci/homebrew-bio

[channels]: https://bioconda.github.io/user/install.html#set-up-channels

[conda]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/

[docker]: https://docs.docker.com/v17.12/install/

[dockerhub]: https://hub.docker.com/r/mbhall88/rasusa

[dpryan79]: https://github.com/dpryan79

[filtlong]: https://github.com/rrwick/Filtlong

[homebrew]: https://docs.brew.sh/Installation

[hyperfine]: https://github.com/sharkdp/hyperfine

[kcov]: https://github.com/SimonKagstrom/kcov

[log-lvl]: https://docs.rs/log/0.4.6/log/enum.Level.html#variants

[mgen-ref]: https://doi.org/10.1099/mgen.0.000294

[pr-help]: https://github.com/bioconda/bioconda-recipes/pull/18690

[pyfastaq]: https://github.com/sanger-pathogens/Fastaq

[quay.io]: https://quay.io/repository/mbhall88/rasusa

[rust]: https://www.rust-lang.org/tools/install

[score]: https://github.com/rrwick/Filtlong#read-scoring

[seed]: https://en.wikipedia.org/wiki/Random_seed

[seqtk]: https://github.com/lh3/seqtk

[singularity]: https://sylabs.io/guides/3.4/user-guide/quick_start.html#quick-installation-steps

[snakemake]: https://snakemake.readthedocs.io/en/stable/

[triples]: https://clang.llvm.org/docs/CrossCompilation.html#target-triple

[wrapper]: https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/rasusa.html

[zstd]: https://github.com/facebook/zstd

