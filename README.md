![rasusa](img/logo.png)

[![Build Status](https://travis-ci.org/mbhall88/rasusa.svg?branch=master)](https://travis-ci.org/mbhall88/rasusa)
[![codecov](https://codecov.io/gh/mbhall88/rasusa/branch/master/graph/badge.svg)](https://codecov.io/gh/mbhall88/rasusa)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![github release version](https://img.shields.io/github/v/release/mbhall88/rasusa)](https://github.com/mbhall88/rasusa/releases)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3546168.svg)](https://doi.org/10.5281/zenodo.3546168)

**Ra**ndomly **su**b**sa**mple sequencing reads to a specified coverage.

[TOC]:#

# Table of Contents
- [Motivation](#motivation)
- [Install](#install)
  - [`cargo`](#cargo)
  - [`conda`](#conda)
  - [`singularity`](#singularity)
  - [`docker`](#docker)
  - [`homebrew`](#homebrew)
  - [Release binaries](#release-binaries)
  - [Build locally](#build-locally)
- [Usage](#usage)
  - [Basic usage](#basic-usage)
  - [Required parameters](#required-parameters)
  - [Optional parameters](#optional-parameters)
  - [Full usage](#full-usage)
  - [Snakemake](#snakemake)
- [Benchmark](#benchmark)
  - [Single long read input](#single-long-read-input)
  - [Paired-end input](#paired-end-input)
- [Contributing](#contributing)
- [Citing](#citing)
  - [Bibtex](#bibtex)


## Motivation

I couldn't find a tool for subsampling reads that met my requirements. All the strategies I could find fell short as they either just wanted a number or  percentage of reads to subsample to or, if they did subsample to a coverage, they assume all reads are the same size (i.e Illumina). As I mostly work with long-read data this posed a problem if I wanted to subsample a file to certain coverage, as length of reads was never taken into account. `rasusa` addresses this shortcoming.

A workaround I had been using for a while was using [`filtlong`][filtlong]. It was simple enough, I just figure out the number of bases I need to achieve a (theoretical) coverage for my sample. Say I have a fastq from an _E. coli_ sample with 5 million reads and I want to subset it to 50x coverage. I just need to multiply the expected size of the sample's genome, 4.6 million base pairs, by the coverage I want and I have my target bases - 230 million base pairs. In `filtlong`, I can do the following

```sh
target=230000000
filtlong --target_bases "$target" reads.fq > reads.50x.fq
```

However, this is technically not the intended function of `filtlong`; it's a quality filtering tool. What you get in the end is a subset of the ["highest scoring"][score] reads at a (theoretical) coverage of 50x. Depending on your circumstances, this might be what you want. However, you bias yourself towards the best/longest reads in the dataset - not a fair representation of your dataset as a whole. There is also the possibility of favouring regions of the genome that produce longer/higher quality reads. [De Maio _et al._][mgen-ref] even found that by randomly subsampling nanopore reads you achieve _better_ genome assemblies than if you had filtered.  

So, depending on your circumstances, an unbiased subsample of your reads might be what you need. And if this is the case, `rasusa` has you covered.

[mgen-ref]: https://doi.org/10.1099/mgen.0.000294

[filtlong]: https://github.com/rrwick/Filtlong

[score]: https://github.com/rrwick/Filtlong#read-scoring

## Install

Some of these installation options require the [`rust` toolchain][rust], which is extremely easy to set up. However, if you do not wish to install `rust` then there are a number of options available.

### `cargo`

[![Crates.io](https://img.shields.io/crates/v/rasusa.svg)](https://crates.io/crates/rasusa)

Prerequisite: [`rust` toolchain][rust]

```sh
cargo install rasusa
```

### `conda`

[![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/rasusa)](https://anaconda.org/bioconda/rasusa)
[![bioconda version](https://anaconda.org/bioconda/rasusa/badges/platforms.svg)](https://anaconda.org/bioconda/rasusa)

Prerequisite: [`conda`][conda] (and bioconda channel [correctly set up][channels])

```sh
conda install rasusa
```

Thank you to Devon Ryan ([@dpryan79][dpryan79]) for [help debugging the bioconda recipe][pr-help].

[channels]: https://bioconda.github.io/user/install.html#set-up-channels
[dpryan79]: https://github.com/dpryan79
[pr-help]: https://github.com/bioconda/bioconda-recipes/pull/18690

### `singularity`

Prerequisite: [`singularity`][singularity]

```sh
URI="docker://mbhall88/rasusa"
singularity exec "$URI" rasusa --help
```

The above will use the latest version. If you want to specify a version then use a [tag][dockerhub] like so.

```sh
VERSION="0.3.0"
URI="docker://mbhall88/rasusa:${VERSION}"
```

### `docker`
[![Docker Image Version (latest semver)](https://img.shields.io/docker/v/mbhall88/rasusa)](https://hub.docker.com/r/mbhall88/rasusa/tags)

Prerequisite: [`docker`][docker]

```sh
docker pull mbhall88/rasusa
docker run mbhall88/rasusa rasusa --help
```

You can find all the available tags on the [Docker Hub repository][dockerhub].

[dockerhub]: https://hub.docker.com/r/mbhall88/rasusa

### `homebrew`

Prerequisite: [`homebrew`][homebrew]

The `homebrew` installation is done via the [homebrew-bio tap][brew-tap].

```sh
brew tap brewsci/bio
brew install rasusa
```

or

```sh
brew install brewsci/bio/rasusa
```

[brew-tap]: https://github.com/brewsci/homebrew-bio

### Release binaries

**tl;dr**: Run the following snippet to download the binary for your system to the current directory and show the help menu.
```shell
OS=$(uname -s)                                                                                                       
if [ "$OS" = "Linux" ]; then                                                                                         
    triple="x86_64-unknown-linux-musl"                                                                              
elif [ "$OS" = "Darwin" ]; then                                                                                        
    triple="x86_64-apple-darwin"                                                         
else                                                      
    echo "ERROR: $OS not a recognised operating system"
fi              
if [ -n "$triple" ]; then   
    URL="https://github.com/mbhall88/rasusa/releases/download/0.3.0/rasusa-0.3.0-${triple}.tar.gz"
    wget "$URL" -O - | tar -xzf -
    ./rasusa --help             
fi
```

These binaries _do not_ require that you have the `rust` toolchain installed.

Currently, there are two pre-compiled binaries available:
- Linux kernel `x86_64-unknown-linux-musl` (works on most Linux distributions I tested)
- OSX kernel `x86_64-apple-darwin` (works for any post-2007 Mac)

An example of downloading one of these binaries using `wget`
```sh
URL="https://github.com/mbhall88/rasusa/releases/download/0.3.0/rasusa-0.3.0-x86_64-unknown-linux-musl.tar.gz"
wget "$URL" -O - | tar -xzf -
./rasusa --help
```

If these binaries do not work on your system please raise an issue and I will potentially
add some additional [target triples][triples].

[triples]: https://clang.llvm.org/docs/CrossCompilation.html#target-triple

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

[rust]: https://www.rust-lang.org/tools/install

[conda]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/

[singularity]: https://sylabs.io/guides/3.4/user-guide/quick_start.html#quick-installation-steps

[docker]: https://docs.docker.com/v17.12/install/

[homebrew]: https://docs.brew.sh/Installation

## Usage

### Basic usage

```
rasusa --input in.fq --coverage 30 --genome-size 4.6mb
```

The above command will output the subsampled file to `stdout`.

Or, if you have paired Illumina

```
rasusa -i r1.fq -i r2.fq --coverage 30 --genome-size 4g -o out.r1.fq -o out.r2.fq
```

For more details on the above options, and additional options, see below.

### Required parameters

There are three required options to run `rasusa`.

#### Input

##### `-i`, `--input`

This option specifies the file(s) containing the reads you would like to subsample.
The file(s) must be valid fasta or fastq format and can be compressed (with a tool such as `gzip`).  
Illumina paired files can be passed in two ways.
1. Using `--input` twice `-i r1.fq -i r2.fq`
2. Using `--input` once, but passing both files immediately after `-i r1.fq r2.fq`

> Bash wizard tip 🧙: Let globs do the work for you `-i r*.fq`

_**Note**: The file format is (lazily) determined from the file name._ File suffixes that are deemed valid are:

-   fasta: `.fa`, `.fasta`, `.fa.gz`, `.fasta.gz`
-   fastq: `.fq`, `.fastq`, `.fq.gz`, `.fastq.gz`

If there is a naming convention you feel is missing, please raise an issue.

#### Coverage

##### `-c`, `--coverage`

This option is used to determine the minimum coverage to subsample the
reads to. It can be specified as an integer (100), a decimal/float
(100.7), or either of the previous suffixed with an 'x' (100x).


_**Note**: Due to the method for determining how many bases are required
to achieve the desired coverage, the actual coverage, in the end, could
be slightly higher than requested. For example, if the last included
read is very long. The log messages should inform you of the actual
coverage in the end._

#### Genome size

##### `-g`, `--genome-size`

The genome size of the input is also required. It is used to determine how many bases are necessary to achieve the desired coverage. This can, of course, be as precise or rough as you like.  
Genome size can be passed in many ways. As a plain old integer (1600), or with a metric suffix (1.6kb). All metric suffixes can have an optional 'b' suffix and be lower, upper, or mixed case. So 'Kb', 'kb' and 'k' would all be inferred as 'kilo'. Valid metric suffixes include:

-   Base (b) - multiplies by 1
-   Kilo (k) - multiplies by 1,000
-   Mega (m) - multiplies by 1,000,000
-   Giga (g) - multiplies by 1,000,000,000
-   Tera (t) - multiplies by 1,000,000,000,000

### Optional parameters

#### Output

##### `-o`, `--output`

NOTE: This parameter is required if passing paired Illumina data.

By default, `rasusa` will output the subsampled file to `stdout` (if one file is given). If you would prefer to specify an output file path, then use this option. If you add the compression suffix `.gz` to the file path then the output will be compressed for you.

Output for Illumina paired files can be specified in the same manner as [`--input`](#input)
1. Using `--output` twice `-o out.r1.fq -o out.r2.fq`
2. Using `--output` once, but passing both files immediately after `-o out.r1.fq out.r2.fq`

The ordering of the output files is assumed to be the same as the input.  
_Note: The output will always be in the same format as the input. You cannot pass fastq as input and ask for fasta as output._

#### Random seed

##### `-s`, `--seed`

This option allows you to specify the [random seed][seed] used by the random subsampler. By explicitly setting this parameter, you make the subsample for the input reproducible. The seed is an integer, and by default it is not set, meaning the operating system will seed the random subsampler. You should only pass this parameter if you are likely to want to subsample the same input file again in the future and want the same subset of reads.

[seed]: https://en.wikipedia.org/wiki/Random_seed

#### Verbosity

##### `-v`

Adding this optional flag will make the logging more verbose. By default, logging will produce messages considered "info" or above (see [here][log-lvl] for more details). If verbosity is switched on, you will additionally get "debug" level logging messages.

[log-lvl]: https://docs.rs/log/0.4.6/log/enum.Level.html#variants

### Full usage

```text
$ rasusa --help

rasusa 0.3.0
Randomly subsample reads to a specified coverage.

USAGE:
    rasusa [FLAGS] [OPTIONS] --coverage <coverage> --genome-size <genome-size> --input <input>...

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information
    -v               Switch on verbosity.

OPTIONS:
    -c, --coverage <coverage>          The desired coverage to sub-sample the reads to.
    -g, --genome-size <genome-size>    Size of the genome to calculate coverage with respect to. i.e 4.3kb, 7Tb, 9000,
                                       4.1MB etc.
    -i, --input <input>...             The fast{a,q} file(s) to subsample. For paired Illumina you may either pass this
                                       flag twice `-i r1.fq -i r2.fq` or give two files consecutively `-i r1.fq r2.fq`.
    -o, --output <output>...           Output file(s), stdout if not present. For paired Illumina you may either pass
                                       this flag twice `-o o1.fq -o o2.fq` or give two files consecutively `-o o1.fq
                                       o2.fq`. NOTE: The order of the pairs is assumed to be the same as that given for
                                       --input. This option is required for paired input.
    -s, --seed <seed>                  Random seed to use.
```

### Snakemake

If you want to use `rasusa` in a [`snakemake`][snakemake] pipeline, it is advised to use the [wrapper][wrapper].

```py
rule subsample:
    input:
        r1="{sample}.r1.fq",
        r2="{sample}.r2.fq",
    output:
        r1="{sample}.subsampled.r1.fq",
        r2="{sample}.subsampled.r2.fq",
    params:
        options="--seed 15",  # optional
        genome_size="3mb",  # required
        coverage=20,  # required
    log:
        "logs/subsample/{sample}.log",
    wrapper:
        "0.70.0/bio/rasusa"
```

*See the [latest wrapper][wrapper] documentation for the most up-to-date version number.*

[snakemake]: https://snakemake.readthedocs.io/en/stable/
[wrapper]: https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/rasusa.html

## Benchmark

> “Time flies like an arrow; fruit flies like a banana.”  
― Anthony G. Oettinger

The real question is: will `rasusa` just needlessly eat away at your precious time on earth?

To do this benchmark, I am going to use [hyperfine][hyperfine].

The data I used comes from

> [Bainomugisa, Arnold, et al. "A complete high-quality MinION nanopore assembly of an extensively drug-resistant Mycobacterium tuberculosis Beijing lineage strain identifies novel variation in repetitive PE/PPE gene regions." Microbial genomics 4.7 (2018).][1]

### Single long read input
Download and rename the fastq
```shell
URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR649/008/SRR6490088/SRR6490088_1.fastq.gz"
wget "$URL" -O - | gzip -d -c > tb.fq
```
The file size is 2.9G, and it has 379,547 reads.  
We benchmark against `filtlong` using the same strategy outlined in [Motivation](#motivation).

```shell
TB_GENOME_SIZE=4411532
COVG=50
TARGET_BASES=$(( TB_GENOME_SIZE * COVG ))
FILTLONG_CMD="filtlong --target_bases $TARGET_BASES tb.fq"
RASUSA_CMD="rasusa -i tb.fq -c $COVG -g $TB_GENOME_SIZE -s 1"
hyperfine --warmup 3 --runs 10 --export-markdown results.md \
     "$FILTLONG_CMD" "$RASUSA_CMD" 
```

#### Results
| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `filtlong --target_bases 220576600 tb.fq` | 111.199 ± 14.926 | 92.716 | 128.849 | 51.64 ± 9.11 |
| `rasusa -i tb.fq -c 50 -g 4411532 -s 1` | 2.153 ± 0.246 | 2.060 | 2.853 | 1.00 |

**Summary**: `rasusa` ran 51.64 ± 9.11 times faster than `filtlong`.


### Paired-end input

Download and then deinterleave the fastq with [`pyfastaq`][pyfastaq]
```shell
URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/008/SRR6488968/SRR6488968.fastq.gz"
wget "$URL" -O - | fastaq deinterleave - r1.fq r2.fq
```
Each file's size is 179M and has 283,590 reads.  
For this benchmark, we will use [`seqtk`][seqtk]. As `seqtk` requires a fixed number of
reads to subsample to, I ran `rasusa` initially and got the number of reads it was using
for its subsample. We will also test `seqtk`'s 2-pass mode as this is analogous to `rasusa`.  
We use a lower coverage for the Illumina as the samples only have ~38x coverage.

```shell
TB_GENOME_SIZE=4411532
COVG=20
NUM_READS=147052
SEQTK_CMD_1="seqtk sample -s 1 r1.fq $NUM_READS > /tmp/r1.fq; seqtk sample -s 1 r2.fq $NUM_READS > /tmp/r2.fq;"
SEQTK_CMD_2="seqtk sample -2 -s 1 r1.fq $NUM_READS > /tmp/r1.fq; seqtk sample -2 -s 1 r2.fq $NUM_READS > /tmp/r2.fq;"
RASUSA_CMD="rasusa -i r1.fq r2.fq -c $COVG -g $TB_GENOME_SIZE -s 1 -o /tmp/r1.fq -o /tmp/r2.fq"
hyperfine --warmup 5 --runs 20 --export-markdown results.md \
     "$SEQTK_CMD_1" "$SEQTK_CMD_2" "$RASUSA_CMD"
```

#### Results

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `seqtk sample -s 1 r1.fq 147052 > /tmp/r1.fq; seqtk sample -s 1 r2.fq 147052 > /tmp/r2.fq;` | 1.435 ± 0.036 | 1.396 | 1.523 | 1.00 |
| `seqtk sample -2 -s 1 r1.fq 147052 > /tmp/r1.fq; seqtk sample -2 -s 1 r2.fq 147052 > /tmp/r2.fq;` | 1.587 ± 0.166 | 1.399 | 2.087 | 1.11 ± 0.12 |
| `rasusa -i r1.fq r2.fq -c 20 -g 4411532 -s 1 -o /tmp/r1.fq -o /tmp/r2.fq` | 1.731 ± 0.079 | 1.611 | 1.982 | 1.21 ± 0.06 |

**Summary**: `seqtk` (1-pass) ran 1.21 ± 0.06 times faster than `rasusa`
and 1.11 ± 0.12 times faster than `seqtk` (2-pass)

So, `rasusa` is just about the same speed as `seqtk` but doesn't require a fixed number of reads - allowing you to avoid doing maths to determine how many reads you need to downsample to a specific coverage. 🤓

[hyperfine]: https://github.com/sharkdp/hyperfine
[1]: https://doi.org/10.1099/mgen.0.000188
[pyfastaq]: https://github.com/sanger-pathogens/Fastaq
[seqtk]: https://github.com/lh3/seqtk

## Contributing

If you would like to help improve `rasusa` you are very welcome!  

For changes to be accepted, they must pass the CI and coverage checks. These include:

-   Code is formatted with `rustfmt`. This can be done by running `cargo fmt` in the project directory.
-   There are no compiler errors/warnings. You can check this by running `cargo clippy --all-features --all-targets -- -D warnings`
-   Code coverage has not reduced. If you want to check coverage before pushing changes, I use [`kcov`][kcov].

[kcov]: https://github.com/SimonKagstrom/kcov

## Citing

If you use `rasusa` in your research, it would be very much appreciated if you could cite it.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3731394.svg)](https://doi.org/10.5281/zenodo.3731394)

> Hall, Michael B. Rasusa: Randomly subsample sequencing reads to a specified coverage. (2019). doi:10.5281/zenodo.3731394

### Bibtex

```Bibtex
@article{
    rasusa2019,
    title={Rasusa: Randomly subsample sequencing reads to a specified coverage},
    DOI={10.5281/zenodo.3731394},
    publisher={Zenodo},
    author={Hall, Michael B.},
    year={2019},
    month={Nov}
}
```
