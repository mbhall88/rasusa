#!/usr/bin/env bash
# A/B benchmark harness for rasusa: wall time + peak RSS per scenario.
#
# Usage:
#   benches/bench.sh run [ref]
#   benches/bench.sh compare <base> [head]
#
#   run [ref]           Build `ref` (a git ref; omit for the current working tree) in
#                        release mode, run every scenario, write
#                        benches/results/<ref>.json.
#   compare <base> [head]
#                        Run `base` and `head` (head defaults to the current working
#                        tree) and print wall-time / peak-RSS ratios. Exits non-zero if
#                        any scenario regressed past the configured thresholds.
#
# Env vars:
#   BENCH_READS         Forwarded to gen_fixtures.sh (default: 1000000)
#   BENCH_SPARSE_READS  Forwarded to gen_fixtures.sh - reads for reads-num-sparse (default: 10000000)
#   BENCH_SPARSE_K      Reads kept (-n) in the reads-num-sparse scenario (default: 1000)
#   BENCH_WARMUP        hyperfine/time-loop warmup runs (default: 3)
#   BENCH_RUNS          hyperfine/time-loop measured runs (default: 10)
#   BENCH_WALL_THRESHOLD  Regression threshold for wall time ratio (default: 1.15)
#   BENCH_RSS_THRESHOLD   Regression threshold for peak RSS ratio (default: 1.15)
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESULTS_DIR="$REPO_ROOT/benches/results"
DATA_DIR="${DATA_DIR:-$REPO_ROOT/data}"
LIB_PY="$REPO_ROOT/benches/lib.py"

WARMUP="${BENCH_WARMUP:-3}"
RUNS="${BENCH_RUNS:-10}"
WALL_THRESHOLD="${BENCH_WALL_THRESHOLD:-1.15}"
RSS_THRESHOLD="${BENCH_RSS_THRESHOLD:-1.15}"
BENCH_SPARSE_READS="${BENCH_SPARSE_READS:-10000000}"
BENCH_SPARSE_K="${BENCH_SPARSE_K:-1000}"
export BENCH_SPARSE_READS

usage() {
    sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'
}

require() {
    command -v "$1" >/dev/null 2>&1 || {
        echo "error: required tool '$1' not found on PATH" >&2
        exit 1
    }
}

sanitize_label() {
    echo "$1" | tr '/' '_'
}

ensure_fixtures() {
    "$REPO_ROOT/benches/gen_fixtures.sh"
    # shellcheck disable=SC1090
    source "$DATA_DIR/meta.env"
}

detect_rss_tool() {
    if [[ "$(uname -s)" == "Darwin" ]]; then
        if [[ -x /usr/bin/time ]]; then
            echo "bsd:/usr/bin/time"
            return
        fi
    fi
    if /usr/bin/time -v true >/dev/null 2>&1; then
        echo "gnu:/usr/bin/time"
        return
    fi
    if command -v gtime >/dev/null 2>&1 && gtime -v true >/dev/null 2>&1; then
        echo "gnu:gtime"
        return
    fi
    echo "none"
}

# Prints peak RSS in bytes for the given command, or nothing if unmeasurable.
measure_peak_rss_bytes() {
    local tool="$1"
    shift
    case "$tool" in
    bsd:*)
        local bin="${tool#bsd:}"
        ("$bin" -l "$@" 2>&1 1>/dev/null) | awk '/maximum resident set size/ {print $1}'
        ;;
    gnu:*)
        local bin="${tool#gnu:}"
        ("$bin" -v "$@" 2>&1 1>/dev/null) | awk -F': ' '/Maximum resident set size/ {print $2 * 1024}'
        ;;
    *)
        echo ""
        ;;
    esac
}

# Builds `ref` in release mode and prints the path to the resulting binary.
# `ref` of "local" (or empty) builds the current working tree in place.
build_ref() {
    local ref="$1"
    if [[ -z "$ref" || "$ref" == "local" ]]; then
        (cd "$REPO_ROOT" && cargo build --release --bin rasusa --quiet) >&2
        echo "$REPO_ROOT/target/release/rasusa"
        return
    fi

    local label
    label="$(sanitize_label "$ref")"
    local wt_dir="$REPO_ROOT/target/bench-worktrees/$label"

    if git -C "$REPO_ROOT" worktree list --porcelain | grep -qx "worktree $wt_dir"; then
        git -C "$wt_dir" checkout --detach --force --quiet "$ref" >&2
    else
        rm -rf "$wt_dir"
        mkdir -p "$(dirname "$wt_dir")"
        git -C "$REPO_ROOT" worktree add --detach --quiet "$wt_dir" "$ref" >&2
    fi

    (cd "$wt_dir" && cargo build --release --bin rasusa --quiet) >&2
    echo "$wt_dir/target/release/rasusa"
}

# Prints a shell-quoted command string for the named scenario.
scenario_cmd_str() {
    local name="$1" bin="$2"
    local n=$((META_BENCH_READS / 4))
    local p=$((META_PAIRED_READS / 2))
    case "$name" in
    reads-num)
        printf '%q ' "$bin" reads "$DATA_DIR/reads.fq" -n "$n" -s 42 -o /dev/null
        ;;
    reads-frac)
        printf '%q ' "$bin" reads "$DATA_DIR/reads.fq" -f 0.25 -s 42 -o /dev/null
        ;;
    reads-coverage)
        printf '%q ' "$bin" reads "$DATA_DIR/reads.fq" -c 30 -g "$META_GENOME_SIZE" -s 42 -o /dev/null
        ;;
    reads-paired)
        printf '%q ' "$bin" reads "$DATA_DIR/r1.fq" "$DATA_DIR/r2.fq" -n "$p" -s 42 -o /dev/null -o /dev/null
        ;;
    reads-num-sparse)
        printf '%q ' "$bin" reads "$DATA_DIR/reads_sparse.fq" -n "$BENCH_SPARSE_K" -s 42 -o /dev/null
        ;;
    aln-stream)
        printf '%q ' "$bin" aln "$REPO_ROOT/tests/cases/test.bam" -c 5 -s 42 --strategy stream -O bam -o /dev/null
        ;;
    aln-fetch)
        printf '%q ' "$bin" aln "$REPO_ROOT/tests/cases/test.bam" -c 5 -s 42 --strategy fetch -O bam -o /dev/null
        ;;
    aln-paired)
        printf '%q ' "$bin" aln "$REPO_ROOT/tests/cases/test.paired.bam" -c 5 -s 42 -O bam -o /dev/null
        ;;
    *)
        echo "unknown scenario: $name" >&2
        return 1
        ;;
    esac
}

scenario_names() {
    echo "reads-num reads-frac reads-coverage reads-paired reads-num-sparse aln-stream aln-fetch aln-paired"
}

run_one_scenario() {
    local name="$1" bin="$2" out_json="$3" rss_tool="$4"
    local cmdstr
    cmdstr="$(scenario_cmd_str "$name" "$bin")"
    local -a cmd
    eval "cmd=( $cmdstr )"

    local mean stddev
    if command -v hyperfine >/dev/null 2>&1; then
        local hjson
        hjson="$(mktemp)"
        hyperfine --warmup "$WARMUP" --runs "$RUNS" --export-json "$hjson" -- "$cmdstr" >&2
        read -r mean stddev < <(python3 "$LIB_PY" parse-hyperfine "$hjson")
        rm -f "$hjson"
    else
        read -r mean stddev < <(python3 "$LIB_PY" time-loop "$WARMUP" "$RUNS" "${cmd[@]}")
    fi

    local rss=""
    if [[ "$rss_tool" != "none" ]]; then
        rss="$(measure_peak_rss_bytes "$rss_tool" "${cmd[@]}")"
    fi

    python3 "$LIB_PY" write-scenario "$out_json" "$name" "$mean" "$stddev" "$rss" >/dev/null
    echo "    wall_mean=${mean}s peak_rss=${rss:-n/a}" >&2
}

cmd_run() {
    local ref="${1:-local}"
    require python3
    require cargo

    ensure_fixtures
    local bin
    bin="$(build_ref "$ref")"

    local label
    label="$(sanitize_label "$ref")"
    mkdir -p "$RESULTS_DIR"
    local out_json="$RESULTS_DIR/$label.json"
    rm -f "$out_json"

    local commit
    if [[ "$ref" == "local" ]]; then
        commit="$(git -C "$REPO_ROOT" rev-parse HEAD)"
        # Tracked-file changes only - untracked scratch files in the working tree
        # shouldn't taint the provenance of an otherwise-clean build.
        if ! git -C "$REPO_ROOT" diff --quiet HEAD --; then
            commit="${commit}+dirty"
        fi
    else
        commit="$(git -C "$REPO_ROOT" rev-parse "$ref")"
    fi

    python3 "$LIB_PY" set-meta "$out_json" ref "$label" >/dev/null
    python3 "$LIB_PY" set-meta "$out_json" commit "$commit" >/dev/null
    python3 "$LIB_PY" set-meta "$out_json" bench_reads "$META_BENCH_READS" >/dev/null

    local rss_tool
    rss_tool="$(detect_rss_tool)"
    if [[ "$rss_tool" == "none" ]]; then
        echo "warning: no peak-RSS measurement tool found (need GNU time -v, gtime, or BSD time -l); RSS will be omitted" >&2
    fi

    for name in $(scenario_names); do
        echo "==> [$label] $name" >&2
        run_one_scenario "$name" "$bin" "$out_json" "$rss_tool"
    done

    echo "$out_json"
}

cmd_compare() {
    local base="${1:?compare requires <base> [head]}"
    local head="${2:-local}"

    local base_json head_json
    base_json="$(cmd_run "$base")"
    head_json="$(cmd_run "$head")"

    python3 "$LIB_PY" compare "$base_json" "$head_json" "$WALL_THRESHOLD" "$RSS_THRESHOLD"
}

main() {
    local sub="${1:-}"
    case "$sub" in
    run)
        shift
        cmd_run "$@"
        ;;
    compare)
        shift
        cmd_compare "$@"
        ;;
    -h | --help | help | "")
        usage
        ;;
    *)
        echo "error: unknown subcommand '$sub'" >&2
        usage >&2
        exit 2
        ;;
    esac
}

main "$@"
