#!/usr/bin/env python3
"""JSON I/O and comparison helpers for benches/bench.sh.

Not meant to be run directly by humans - it's a small dispatch table of subcommands
bench.sh shells out to, so bench.sh doesn't need to hand-roll JSON parsing/formatting
in bash.
"""
import json
import subprocess
import sys
import time


def _load(path):
    try:
        with open(path) as f:
            return json.load(f)
    except FileNotFoundError:
        return {"scenarios": {}}


def _dump(path, data):
    with open(path, "w") as f:
        json.dump(data, f, indent=2, sort_keys=True)
        f.write("\n")


def cmd_write_scenario(args):
    # write-scenario <result.json> <scenario> <wall_mean_s> <wall_stddev_s> <peak_rss_bytes|"">
    path, scenario, wall_mean, wall_stddev, rss = args
    data = _load(path)
    data.setdefault("scenarios", {})[scenario] = {
        "wall_mean_s": float(wall_mean),
        "wall_stddev_s": float(wall_stddev),
        "peak_rss_bytes": int(float(rss)) if rss not in ("", "null") else None,
    }
    _dump(path, data)


def cmd_set_meta(args):
    # set-meta <result.json> <key> <value>
    path, key, value = args
    data = _load(path)
    data[key] = value
    _dump(path, data)


def cmd_parse_hyperfine(args):
    # parse-hyperfine <hyperfine-export.json>  -> prints "<mean> <stddev>"
    (path,) = args
    with open(path) as f:
        data = json.load(f)
    r = data["results"][0]
    print(r["mean"], r["stddev"])


def cmd_time_loop(args):
    # time-loop <warmup> <runs> <cmd...>  -> prints "<mean> <stddev>"
    warmup = int(args[0])
    runs = int(args[1])
    cmd = args[2:]
    for _ in range(warmup):
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    times = []
    for _ in range(runs):
        t0 = time.perf_counter()
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        times.append(time.perf_counter() - t0)
    mean = sum(times) / len(times)
    variance = sum((t - mean) ** 2 for t in times) / len(times)
    print(mean, variance**0.5)


def cmd_compare(args):
    # compare <base.json> <head.json> <wall_threshold> <rss_threshold>
    base_path, head_path, wall_threshold, rss_threshold = args
    wall_threshold = float(wall_threshold)
    rss_threshold = float(rss_threshold)
    base = _load(base_path)
    head = _load(head_path)

    base_scn = base.get("scenarios", {})
    head_scn = head.get("scenarios", {})
    names = sorted(set(base_scn) | set(head_scn))

    header = (
        "scenario",
        "base wall(s)",
        "head wall(s)",
        "wall ratio",
        "base RSS(MB)",
        "head RSS(MB)",
        "RSS ratio",
    )
    rows = [header]
    regressed = []

    for name in names:
        b = base_scn.get(name)
        h = head_scn.get(name)
        if not b or not h:
            rows.append((name, "missing", "-", "-", "-", "-", "-"))
            continue

        wall_ratio = h["wall_mean_s"] / b["wall_mean_s"] if b["wall_mean_s"] else float("nan")

        rss_ratio = None
        if b.get("peak_rss_bytes") and h.get("peak_rss_bytes"):
            rss_ratio = h["peak_rss_bytes"] / b["peak_rss_bytes"]

        rows.append((
            name,
            f"{b['wall_mean_s']:.3f}",
            f"{h['wall_mean_s']:.3f}",
            f"{wall_ratio:.2f}",
            f"{b['peak_rss_bytes'] / 1e6:.1f}" if b.get("peak_rss_bytes") else "-",
            f"{h['peak_rss_bytes'] / 1e6:.1f}" if h.get("peak_rss_bytes") else "-",
            f"{rss_ratio:.2f}" if rss_ratio is not None else "-",
        ))

        if wall_ratio > wall_threshold:
            regressed.append(f"{name}: wall time {wall_ratio:.2f}x base (threshold {wall_threshold}x)")
        if rss_ratio is not None and rss_ratio > rss_threshold:
            regressed.append(f"{name}: peak RSS {rss_ratio:.2f}x base (threshold {rss_threshold}x)")

    widths = [max(len(str(r[i])) for r in rows) for i in range(len(header))]

    def fmt_row(r):
        return "  ".join(str(v).ljust(widths[i]) for i, v in enumerate(r))

    print(fmt_row(rows[0]))
    print("  ".join("-" * w for w in widths))
    for r in rows[1:]:
        print(fmt_row(r))

    if regressed:
        print("\nREGRESSIONS DETECTED:", file=sys.stderr)
        for line in regressed:
            print(f"  - {line}", file=sys.stderr)
        sys.exit(1)


COMMANDS = {
    "write-scenario": cmd_write_scenario,
    "set-meta": cmd_set_meta,
    "parse-hyperfine": cmd_parse_hyperfine,
    "time-loop": cmd_time_loop,
    "compare": cmd_compare,
}


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        print(f"usage: {sys.argv[0]} <{'|'.join(COMMANDS)}> [args...]", file=sys.stderr)
        sys.exit(2)
    COMMANDS[sys.argv[1]](sys.argv[2:])


if __name__ == "__main__":
    main()
