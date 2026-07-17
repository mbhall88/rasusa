#!/usr/bin/env bash
# Deterministically generate synthetic FASTQ fixtures for the benchmark harness.
#
# Usage: benches/gen_fixtures.sh
#
# Env vars:
#   BENCH_READS  - number of single-end reads to generate (default: 1000000)
#   READ_LEN     - length (bp) of each generated read (default: 150)
#   DATA_DIR     - output directory (default: <repo>/data)
#
# Output is a deterministic function of BENCH_READS and READ_LEN alone (fixed internal
# seed), so regenerating with the same values produces byte-identical files. Files are
# cached and only regenerated if those values change or files are missing.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BENCH_READS="${BENCH_READS:-1000000}"
READ_LEN="${READ_LEN:-150}"
DATA_DIR="${DATA_DIR:-$REPO_ROOT/data}"

READS_FQ="$DATA_DIR/reads.fq"
R1_FQ="$DATA_DIR/r1.fq"
R2_FQ="$DATA_DIR/r2.fq"
META="$DATA_DIR/meta.env"

mkdir -p "$DATA_DIR"

if [[ -f "$META" ]]; then
    # shellcheck disable=SC1090
    source "$META"
    if [[ "${META_BENCH_READS:-}" == "$BENCH_READS" && "${META_READ_LEN:-}" == "$READ_LEN" \
        && -f "$READS_FQ" && -f "$R1_FQ" && -f "$R2_FQ" ]]; then
        echo "[gen_fixtures] cached fixtures for BENCH_READS=$BENCH_READS READ_LEN=$READ_LEN already exist in $DATA_DIR" >&2
        exit 0
    fi
fi

echo "[gen_fixtures] generating fixtures for BENCH_READS=$BENCH_READS READ_LEN=$READ_LEN into $DATA_DIR" >&2

PAIRED_READS=$(( BENCH_READS / 5 ))
if [[ "$PAIRED_READS" -lt 2 ]]; then
    PAIRED_READS=2
fi

# Both the single-end and paired-end sets are generated in one awk process, sharing the
# same PRNG/sequence-generation functions (a Park-Miller LCG seeded with fixed
# constants, so output only depends on BENCH_READS/READ_LEN/PAIRED_READS, never on the
# host's random source or awk implementation's built-in rand()). Each set reseeds `state`
# independently so they don't affect each other.
awk -v n="$BENCH_READS" -v pn="$PAIRED_READS" -v len="$READ_LEN" \
    -v seed1=42 -v seed2=142857 \
    -v out="$READS_FQ" -v out1="$R1_FQ" -v out2="$R2_FQ" '
function rnd() {
    state = (state * 48271) % 2147483647
    return state / 2147483647
}
function randseq(   seq, j, idx) {
    seq = ""
    for (j = 0; j < len; j++) {
        idx = int(rnd() * 4) + 1
        seq = seq substr(bases, idx, 1)
    }
    return seq
}
BEGIN {
    bases = "ACGT"
    qual = sprintf("%*s", len, "")
    gsub(/ /, "I", qual)

    state = seed1
    for (i = 0; i < n; i++) {
        printf "@read%d\n%s\n+\n%s\n", i, randseq(), qual > out
    }

    state = seed2
    for (i = 0; i < pn; i++) {
        printf "@pair%d/1\n%s\n+\n%s\n", i, randseq(), qual > out1
        printf "@pair%d/2\n%s\n+\n%s\n", i, randseq(), qual > out2
    }
}'

TOTAL_BASES=$(( BENCH_READS * READ_LEN ))
# Pick a genome size such that a 30x coverage target uses well under the total bases
# generated, so the `reads-coverage` scenario is always satisfiable.
GENOME_SIZE=$(( TOTAL_BASES / 60 ))

cat > "$META" << EOF
META_BENCH_READS=$BENCH_READS
META_READ_LEN=$READ_LEN
META_TOTAL_BASES=$TOTAL_BASES
META_GENOME_SIZE=$GENOME_SIZE
META_PAIRED_READS=$PAIRED_READS
EOF

echo "[gen_fixtures] done: $BENCH_READS single-end reads ($TOTAL_BASES bp), $PAIRED_READS read pairs" >&2
