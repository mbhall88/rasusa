#!/usr/bin/env python3
"""Helpers for benches/update_readme.sh: render benchmark blocks and splice them into
README.md between the `<!-- BENCH:<name>:START/END -->` markers.
"""
import json
import re
import sys


def _relative_str(mean, stddev, fastest_mean, fastest_stddev, with_stddev=True):
    ratio = mean / fastest_mean
    if with_stddev:
        rel_err = ratio * ((stddev / mean) ** 2 + (fastest_stddev / fastest_mean) ** 2) ** 0.5
        return f"{ratio:.2f} ± {rel_err:.2f}"
    return f"{ratio:.2f}"


def _load_results(json_path):
    with open(json_path) as f:
        return json.load(f)["results"]


def _read_table(md_path):
    with open(md_path) as f:
        return f.read().rstrip("\n")


def cmd_single_block(args):
    md_path, json_path = args
    results = _load_results(json_path)
    rasusa = next(r for r in results if "rasusa" in r["command"])
    other = next(r for r in results if r is not rasusa)

    table = _read_table(md_path)
    relative = _relative_str(other["mean"], other["stddev"], rasusa["mean"], rasusa["stddev"])
    print(table)
    print()
    print(f"**Summary**: `rasusa` ran {relative} times faster than `filtlong`.")


def cmd_paired_block(args):
    md_path, json_path = args
    results = _load_results(json_path)
    rasusa = next(r for r in results if "rasusa" in r["command"])
    others = [r for r in results if r is not rasusa]
    # Order: 1-pass seqtk, 2-pass seqtk, matching the hyperfine invocation order.
    ratio1 = _relative_str(others[0]["mean"], others[0]["stddev"], rasusa["mean"], rasusa["stddev"], with_stddev=False)
    ratio2 = _relative_str(others[1]["mean"], others[1]["stddev"], rasusa["mean"], rasusa["stddev"], with_stddev=False)

    table = _read_table(md_path)
    print(table)
    print()
    print(
        f"**Summary**: `rasusa reads` ran {ratio1} times faster than `seqtk` (1-pass) "
        f"and {ratio2} times faster than `seqtk` (2-pass)"
    )


def cmd_splice(args):
    readme_path, marker, content_path = args
    with open(readme_path) as f:
        readme = f.read()
    with open(content_path) as f:
        content = f.read().rstrip("\n")

    start = f"<!-- BENCH:{marker}:START -->"
    end = f"<!-- BENCH:{marker}:END -->"
    pattern = re.compile(re.escape(start) + r".*?" + re.escape(end), re.DOTALL)
    if not pattern.search(readme):
        print(f"error: markers {start!r}/{end!r} not found in {readme_path}", file=sys.stderr)
        sys.exit(1)

    replacement = f"{start}\n{content}\n{end}"
    readme = pattern.sub(lambda _: replacement, readme, count=1)
    with open(readme_path, "w") as f:
        f.write(readme)


COMMANDS = {
    "single-block": cmd_single_block,
    "paired-block": cmd_paired_block,
    "splice": cmd_splice,
}


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:
        print(f"usage: {sys.argv[0]} <{'|'.join(COMMANDS)}> [args...]", file=sys.stderr)
        sys.exit(2)
    COMMANDS[sys.argv[1]](sys.argv[2:])


if __name__ == "__main__":
    main()
