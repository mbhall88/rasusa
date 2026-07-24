#!/usr/bin/env bash
# Regenerates the top-level README's tool-comparison benchmark tables (rasusa vs
# filtlong/seqtk) and splices them between the `<!-- BENCH:*:START/END -->` markers.
#
# Downloads two real public sequencing datasets on first run (~2.9 GB total) and caches
# them under data/download-cache/ (override with BENCH_DOWNLOAD_CACHE). Intended to run
# from .github/workflows/benchmark.yaml on release; not routine local use.
#
# Requires on PATH: hyperfine, filtlong, seqtk, wget, python3, cargo (unless
# RASUSA_BIN is set to a prebuilt binary). filtlong/seqtk are available via bioconda:
# `conda install -c bioconda filtlong seqtk`. The paired-end dataset is interleaved -
# deinterleaving it is a one-off data-prep step, done with a plain awk split below
# rather than pulling in a tool dependency for it (see the comment at that step for why).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
README="$REPO_ROOT/README.md"
CACHE_DIR="${BENCH_DOWNLOAD_CACHE:-$REPO_ROOT/data/download-cache}"
WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

mkdir -p "$CACHE_DIR"

require() {
    command -v "$1" >/dev/null 2>&1 || {
        echo "error: required tool '$1' not found on PATH" >&2
        exit 1
    }
}
require hyperfine
require filtlong
require seqtk
require wget
require python3

if [[ -z "${RASUSA_BIN:-}" ]]; then
    echo "[update_readme] building rasusa --release..." >&2
    require cargo
    (cd "$REPO_ROOT" && cargo build --release --bin rasusa --quiet)
    RASUSA_BIN="$REPO_ROOT/target/release/rasusa"
fi
[[ -x "$RASUSA_BIN" ]] || {
    echo "error: RASUSA_BIN '$RASUSA_BIN' is not an executable" >&2
    exit 1
}

# Run benchmarked commands from CACHE_DIR and reference input FASTQs by basename, so the
# rendered README table shows readable commands (`rasusa reads tb.fq ...`) instead of the
# CI runner's full checkout path.
cd "$CACHE_DIR"

### Single long-read input (rasusa vs filtlong) ###

TB_FQ="$CACHE_DIR/tb.fq"
if [[ ! -f "$TB_FQ" ]]; then
    echo "[update_readme] downloading single long-read dataset..." >&2
    URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR649/008/SRR6490088/SRR6490088_1.fastq.gz"
    wget -q "$URL" -O - | gzip -d -c > "$TB_FQ"
else
    echo "[update_readme] using cached $TB_FQ" >&2
fi

TB_GENOME_SIZE=4411532
COVG=50
TARGET_BASES=$((TB_GENOME_SIZE * COVG))
TB_FQ_BASENAME="$(basename "$TB_FQ")"
FILTLONG_CMD="filtlong --target_bases $TARGET_BASES $TB_FQ_BASENAME"
RASUSA_CMD="$RASUSA_BIN reads $TB_FQ_BASENAME -c $COVG -g $TB_GENOME_SIZE -s 1 -o /dev/null"

echo "[update_readme] running single-input benchmark..." >&2
hyperfine --warmup 3 --runs 10 \
    --export-markdown "$WORK_DIR/single.md" --export-json "$WORK_DIR/single.json" \
    "$FILTLONG_CMD" "$RASUSA_CMD"

python3 "$REPO_ROOT/benches/render_readme.py" single-block "$WORK_DIR/single.md" "$WORK_DIR/single.json" > "$WORK_DIR/single-block.md"
python3 "$REPO_ROOT/benches/render_readme.py" splice "$README" single "$WORK_DIR/single-block.md"

### Paired-end input (rasusa vs seqtk) ###

R1_FQ="$CACHE_DIR/r1.fq"
R2_FQ="$CACHE_DIR/r2.fq"
if [[ ! -f "$R1_FQ" || ! -f "$R2_FQ" ]]; then
    echo "[update_readme] downloading paired-end dataset..." >&2
    URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR648/008/SRR6488968/SRR6488968.fastq.gz"
    # Deinterleave (alternating 4-line R1/R2 records) with awk rather than a tool
    # dependency (e.g. pyfastaq's `fastaq deinterleave`) - this is a one-off data-prep
    # step, not part of what's being benchmarked, so it doesn't need to match any
    # particular tool's implementation.
    wget -q "$URL" -O - | gzip -d -c - | awk -v r1="$R1_FQ" -v r2="$R2_FQ" '{
        if (((NR - 1) % 8) < 4) print > r1
        else print > r2
    }'
else
    echo "[update_readme] using cached $R1_FQ / $R2_FQ" >&2
fi

NUM_READS=140000
R1_FQ_BASENAME="$(basename "$R1_FQ")"
R2_FQ_BASENAME="$(basename "$R2_FQ")"
SEQTK_CMD_1="seqtk sample -s 1 $R1_FQ_BASENAME $NUM_READS > $WORK_DIR/o1.fq; seqtk sample -s 1 $R2_FQ_BASENAME $NUM_READS > $WORK_DIR/o2.fq;"
SEQTK_CMD_2="seqtk sample -2 -s 1 $R1_FQ_BASENAME $NUM_READS > $WORK_DIR/o1.fq; seqtk sample -2 -s 1 $R2_FQ_BASENAME $NUM_READS > $WORK_DIR/o2.fq;"
RASUSA_CMD="$RASUSA_BIN reads $R1_FQ_BASENAME $R2_FQ_BASENAME -n $NUM_READS -s 1 -o $WORK_DIR/o1.fq -o $WORK_DIR/o2.fq"

echo "[update_readme] running paired-input benchmark..." >&2
hyperfine --warmup 10 --runs 100 \
    --export-markdown "$WORK_DIR/paired.md" --export-json "$WORK_DIR/paired.json" \
    "$SEQTK_CMD_1" "$SEQTK_CMD_2" "$RASUSA_CMD"

python3 "$REPO_ROOT/benches/render_readme.py" paired-block "$WORK_DIR/paired.md" "$WORK_DIR/paired.json" > "$WORK_DIR/paired-block.md"
python3 "$REPO_ROOT/benches/render_readme.py" splice "$README" paired "$WORK_DIR/paired-block.md"

echo "[update_readme] README.md benchmark tables updated" >&2
