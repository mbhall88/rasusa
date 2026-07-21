# Benchmarks

Two independent things live here:

1. **`bench.sh`** - an internal A/B regression harness. Compares two git refs (or the
   current working tree) on wall time and peak RSS, for a fixed set of scenarios. Run it
   locally to prove a change didn't regress performance before opening a PR. It also runs
   automatically on every PR via `.github/workflows/benchmark-pr.yaml`
   (`compare origin/main`), which posts the comparison table to the job summary. That
   check is intentionally non-blocking (advisory only) - a noisy shared CI runner isn't a
   reliable enough signal to gate merges on.
2. **`update_readme.sh`** - regenerates the tool-comparison tables in the top-level
   `README.md` (`rasusa` vs `filtlong`/`seqtk`). This downloads real public datasets and
   is run from a release-triggered GitHub Actions workflow
   (`.github/workflows/benchmark.yaml`), not typically by hand.

## `bench.sh`

### Prerequisites

- [`hyperfine`](https://github.com/sharkdp/hyperfine) for wall-time measurement.
  Install via `brew install hyperfine` or `cargo install hyperfine`. If absent, the
  script falls back to a portable Python timing loop with the same warmup/runs
  semantics - slightly less precise, but no hard dependency.
- Peak RSS is measured with whatever is available:
  - macOS: BSD `/usr/bin/time -l` (built in).
  - Linux: GNU `time -v` (usually `/usr/bin/time`, package `time`), or
    [`gtime`](https://formulae.brew.sh/formula/gtime) (`brew install gnu-time`) if
    running the Linux-style tool via Homebrew on macOS.
  - If none are found, RSS is omitted (wall time is still measured) and a warning is
    printed.
- `python3` for JSON handling (`benches/lib.py`) - no third-party packages required.
- `cargo`/`rustc`, obviously.

### Usage

```shell
# Benchmark the current working tree (uncommitted changes included):
benches/bench.sh run

# Benchmark a specific git ref (branch, tag, or commit):
benches/bench.sh run origin/main

# Compare two refs head-to-head and fail if either regressed:
benches/bench.sh compare origin/main          # base=origin/main, head=current working tree
benches/bench.sh compare origin/main HEAD~1   # base=origin/main, head=HEAD~1
```

`run [ref]` builds `ref` in `--release` mode (a git ref is built in an isolated
`git worktree` under `target/bench-worktrees/`; omitting `ref`, or passing `local`,
builds the current working tree in place), runs every scenario, and writes
`benches/results/<ref>.json` (`ref` is slash-sanitised for the filename). It prints the
path to the written JSON file on stdout.

`compare <base> [head]` runs `run` for both refs and prints a table of wall-time and
peak-RSS ratios (`head / base`). It exits non-zero if any scenario's ratio exceeds the
configured threshold - i.e. a regression.

`benches/results/` is gitignored - it's scratch output, not a reference. The one
committed exception is **`benches/baseline.json`**, a snapshot of `bench.sh run` on
`main` at a point in time, kept for reference. It is *not* used as an automatic CI gate;
use `compare` against a real base ref (e.g. `origin/main`) for that.

### Fixtures

`benches/gen_fixtures.sh` deterministically generates synthetic FASTQ into `data/`
(gitignored) - a single-end file and a paired-end pair, uniform-random ACGT bases from a
fixed-seed PRNG, so the same `BENCH_READS` always produces byte-identical files. It's
called automatically by `bench.sh run`/`compare` and caches its output
(`data/meta.env` records the parameters used), so it only regenerates when
`BENCH_READS`/`READ_LEN` change.

```shell
BENCH_READS=1000000 benches/gen_fixtures.sh   # default: 1,000,000 single-end reads
```

Alignment scenarios use the small, already-committed BAM fixtures under
`tests/cases/` (`test.bam`, `test.paired.bam`) rather than generated data.

### Scenarios

| Scenario         | Command                                                         |
|-------------------|------------------------------------------------------------------|
| `reads-num`       | `reads` on generated single-end FASTQ, `-n <count>`               |
| `reads-frac`      | `reads` on generated single-end FASTQ, `-f 0.25`                  |
| `reads-coverage`  | `reads` on generated single-end FASTQ, `-c 30 -g <genome size>`   |
| `reads-paired`    | `reads` on generated paired-end FASTQ, `-n <count>`                |
| `aln-stream`      | `aln` on `tests/cases/test.bam`, `--strategy stream`               |
| `aln-fetch`       | `aln` on `tests/cases/test.bam`, `--strategy fetch`                 |
| `aln-paired`      | `aln` on `tests/cases/test.paired.bam` (paired-end alignment)       |

All scenarios use a fixed `--seed 42` (or `142857` internally for paired fixture
generation) so results are reproducible run-to-run modulo real timing/memory noise.

### Env vars

| Var                     | Default   | Meaning                                      |
|--------------------------|-----------|-----------------------------------------------|
| `BENCH_READS`            | `1000000` | Reads generated for the `reads-*` scenarios   |
| `BENCH_WARMUP`           | `3`       | Warmup runs before measurement                |
| `BENCH_RUNS`             | `10`      | Measured runs                                 |
| `BENCH_WALL_THRESHOLD`   | `1.15`    | `compare` fails if wall-time ratio exceeds this |
| `BENCH_RSS_THRESHOLD`    | `1.15`    | `compare` fails if peak-RSS ratio exceeds this  |

### Result JSON schema

```jsonc
{
  "ref": "origin_main",
  "commit": "<sha>[+dirty]",
  "bench_reads": "1000000",
  "scenarios": {
    "<scenario-name>": {
      "wall_mean_s": 0.159,
      "wall_stddev_s": 0.002,
      "peak_rss_bytes": 13058048
    },
    "...": { }
  }
}
```

## `update_readme.sh`

See the top-level `README.md` [Benchmark](../README.md#benchmark) section for what this
regenerates and why. Run via the `benchmark.yaml` workflow on release (or
`workflow_dispatch`); not intended for routine local use since it downloads ~2.9 GB of
public sequencing data (cached between runs where possible).
