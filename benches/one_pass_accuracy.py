#!/usr/bin/env python3
"""Measures how closely `rasusa reads --one-pass -f <fraction>` sticks to the
fraction you asked for, across a range of input sizes and fractions.

## Why this exists, in plain terms

`--one-pass` decides whether to keep each read with a coin flip weighted by the
requested fraction, made independently for every read, in a single streaming pass.
That's what makes it fast (see the `reads-frac` vs `reads-frac-one-pass` scenarios in
`bench.sh`) - but it also means the number of reads it actually keeps is not exactly
`fraction * total_reads`, the way the default two-pass mode's *is*. It's a random
outcome that lands close to the target on average, with some spread around it - the
same way flipping a coin 100 times doesn't always land exactly 50 heads.

That spread shrinks as the input gets bigger: quadruple the number of reads and the
typical spread roughly halves. This script quantifies that by actually running the
`rasusa` binary many times at different input sizes and fractions, and comparing what
it kept against what was asked for.

## Usage

    cargo build --release --bin rasusa   # if not already built
    benches/one_pass_accuracy.py [--sizes N,N,...] [--fractions P,P,...] [--reps N]

Prints progress to stderr and a markdown results table to stdout. Takes roughly a
minute at the defaults (5 sizes x 8 fractions x 30 reps = 1200 runs).
"""
import argparse
import itertools
import re
import statistics
import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
BIN = REPO_ROOT / "target" / "release" / "rasusa"
SOURCE_FQ = REPO_ROOT / "data" / "reads.fq"
GEN_FIXTURES = REPO_ROOT / "benches" / "gen_fixtures.sh"

DEFAULT_SIZES = [100, 1_000, 10_000, 100_000, 1_000_000]
DEFAULT_FRACTIONS = [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99]
DEFAULT_REPS = 30

REALISED_RE = re.compile(r"realised fraction: ([\d.]+)\)")


def ensure_fixture(min_reads):
    """Makes sure data/reads.fq exists and has at least `min_reads` records,
    (re)generating it via gen_fixtures.sh if not."""
    have = 0
    if SOURCE_FQ.exists():
        with open(SOURCE_FQ, "rb") as f:
            have = sum(1 for _ in f) // 4
    if have >= min_reads:
        return
    subprocess.run(
        [str(GEN_FIXTURES)],
        env={**os.environ, "BENCH_READS": str(min_reads)},
        check=True,
    )


def write_subset(size, tmp_dir):
    """Writes the first `size` records of data/reads.fq to a fixture in `tmp_dir`."""
    subset_path = Path(tmp_dir) / f"reads_{size}.fq"
    with open(SOURCE_FQ) as src, open(subset_path, "w") as dst:
        for line in itertools.islice(src, size * 4):
            dst.write(line)
    return subset_path


def realised_fraction(fixture, fraction, seed):
    proc = subprocess.run(
        [str(BIN), "reads", str(fixture), "-f", str(fraction), "--one-pass",
         "-s", str(seed), "-v", "-o", "/dev/null"],
        capture_output=True, text=True, check=True,
    )
    match = REALISED_RE.search(proc.stderr)
    if not match:
        raise RuntimeError(f"no realised fraction in stderr: {proc.stderr!r}")
    return float(match.group(1))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sizes", default=",".join(str(n) for n in DEFAULT_SIZES),
                         help="comma-separated input sizes (read counts)")
    parser.add_argument("--fractions", default=",".join(str(p) for p in DEFAULT_FRACTIONS),
                         help="comma-separated requested fractions")
    parser.add_argument("--reps", type=int, default=DEFAULT_REPS,
                         help="runs per (size, fraction) combination, each with a different seed")
    args = parser.parse_args()

    sizes = [int(s) for s in args.sizes.split(",")]
    fractions = [float(p) for p in args.fractions.split(",")]

    if not BIN.exists():
        sys.exit(f"error: {BIN} not found - run `cargo build --release --bin rasusa` first")

    ensure_fixture(max(sizes))

    rows = []
    with tempfile.TemporaryDirectory() as tmp_dir:
        for size in sizes:
            fixture = write_subset(size, tmp_dir)
            for fraction in fractions:
                realised = [realised_fraction(fixture, fraction, seed) for seed in range(args.reps)]
                mean = statistics.mean(realised)
                empirical_std = statistics.pstdev(realised)
                max_abs_dev = max(abs(r - fraction) for r in realised)
                theoretical_std = (fraction * (1 - fraction) / size) ** 0.5
                rows.append((size, fraction, mean, empirical_std, theoretical_std, max_abs_dev))
                print(
                    f"n={size:>8} p={fraction:<5} mean={mean:.4f} std={empirical_std:.5f} "
                    f"theory_std={theoretical_std:.5f} max|dev|={max_abs_dev:.5f}",
                    file=sys.stderr,
                )

    print("\n| n | requested fraction | mean realised fraction | observed spread (std) | "
          "expected spread (std) | worst deviation seen |")
    print("|---|---|---|---|---|---|")
    for size, fraction, mean, empirical_std, theoretical_std, max_abs_dev in rows:
        print(f"| {size:,} | {fraction} | {mean:.4f} | {empirical_std:.5f} | "
              f"{theoretical_std:.5f} | {max_abs_dev:.5f} |")


if __name__ == "__main__":
    main()
