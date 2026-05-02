#!/usr/bin/env python3
"""Diff per-model wall-time between two feature_coverage_report.md files.

Reads the markdown report emitted by benchmark_feature_coverage.py — one
from a baseline (e.g. main) and one from a candidate (e.g. PR HEAD) —
and prints a sorted slowdown table.

Intended consumer: the perf-diff CI job, which runs both branches on the
same runner and uploads the output as a PR artifact.  Not a hard gate:
shared CI runners are noisy enough that single-model deltas of 30%+
happen from neighbour-VM contention, not real regressions.  The diff is
an aid to human review, not a blocker.

Output columns:
  model                  model name
  base_s                 wall-time on baseline (seconds, last-rep total)
  head_s                 wall-time on candidate
  delta_s                head_s - base_s
  delta_pct              (delta_s / base_s) * 100, signed
  flag                   blank / SLOWER / FASTER / NEW / GONE

Sorting: largest absolute delta_pct first.

Usage:
  perf_diff.py <base_report.md> <head_report.md> [--threshold 15]

The --threshold argument controls the SLOWER / FASTER tag boundary
(default 15%).  Models inside the threshold get a blank flag.  Exit
status is always 0; the script reports findings, it does not gate.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

# `### model_name` followed within a few lines by:
#   `- RM reps: N, wall time: X.XXXs`
_MODEL_HEADER = re.compile(r"^###\s+(\S+)")
_WALL_LINE = re.compile(r"^-\s+RM reps:.*wall time:\s+([0-9.]+)s")


def parse_report(path: Path) -> dict[str, float]:
    """Return {model_name: wall_time_seconds} from a benchmark report."""
    out: dict[str, float] = {}
    current: str | None = None
    for line in path.read_text().splitlines():
        m = _MODEL_HEADER.match(line)
        if m:
            current = m.group(1)
            continue
        if current is None:
            continue
        w = _WALL_LINE.match(line)
        if w:
            out[current] = float(w.group(1))
            current = None  # consume; ignore later wall-time lines
    return out


def diff(base: dict[str, float], head: dict[str, float], threshold_pct: float):
    rows = []
    all_models = sorted(set(base) | set(head))
    for m in all_models:
        b = base.get(m)
        h = head.get(m)
        if b is None and h is not None:
            rows.append((m, None, h, None, None, "NEW"))
        elif b is not None and h is None:
            rows.append((m, b, None, None, None, "GONE"))
        else:
            assert b is not None and h is not None
            d = h - b
            pct = (d / b * 100.0) if b > 0 else 0.0
            if pct > threshold_pct:
                flag = "SLOWER"
            elif pct < -threshold_pct:
                flag = "FASTER"
            else:
                flag = ""
            rows.append((m, b, h, d, pct, flag))
    # Sort by abs(pct) desc; missing values sink to bottom.
    rows.sort(key=lambda r: abs(r[4]) if r[4] is not None else -1, reverse=True)
    return rows


def render(rows, threshold_pct: float) -> str:
    out = []
    out.append(f"Per-model wall-time delta (threshold ±{threshold_pct:.0f}%)")
    out.append("=" * 72)
    out.append(f"{'model':<40} {'base_s':>9} {'head_s':>9} {'Δs':>8} {'Δ%':>7}  flag")
    out.append("-" * 72)
    n_slower = n_faster = 0
    for m, b, h, d, pct, flag in rows:
        bs = f"{b:.3f}" if b is not None else "  -"
        hs = f"{h:.3f}" if h is not None else "  -"
        ds = f"{d:+.3f}" if d is not None else "    -"
        ps = f"{pct:+.1f}" if pct is not None else "    -"
        out.append(f"{m:<40} {bs:>9} {hs:>9} {ds:>8} {ps:>7}  {flag}")
        if flag == "SLOWER":
            n_slower += 1
        elif flag == "FASTER":
            n_faster += 1
    out.append("-" * 72)
    out.append(
        f"{n_slower} model(s) >{threshold_pct:.0f}% slower, "
        f"{n_faster} model(s) >{threshold_pct:.0f}% faster"
    )
    out.append("")
    out.append("Note: shared-runner wall-time noise is typically ±20-30% per model.")
    out.append("Treat single-model deltas with skepticism; correlated movement")
    out.append("across many models is the real signal.")
    return "\n".join(out)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("base_report", type=Path)
    p.add_argument("head_report", type=Path)
    p.add_argument(
        "--threshold",
        type=float,
        default=15.0,
        help="percent change to flag as SLOWER/FASTER (default 15)",
    )
    args = p.parse_args()

    base = parse_report(args.base_report)
    head = parse_report(args.head_report)
    if not base:
        print(f"warning: no models parsed from {args.base_report}", file=sys.stderr)
    if not head:
        print(f"warning: no models parsed from {args.head_report}", file=sys.stderr)
    print(render(diff(base, head, args.threshold), args.threshold))
    return 0


if __name__ == "__main__":
    sys.exit(main())
