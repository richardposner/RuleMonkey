#!/usr/bin/env python3
"""basicmodels harness — RM vs NFsim 100-rep ensemble parity check.

Thin wrapper around `harness/benchmark_full.py` that points its verdict
machinery at `tests/reference/basicmodels/` instead of the corpus
references. Same z-score / `tz_max < T_model` math; same report shape.

The 29 BNGL files under `tests/models/nfsim_basicmodels/` are a parity
test set derived from NFsim's `validate/basicModels/` regression suite.
Most models are small and fast; the full suite runs in a few minutes.
See `tests/reference/basicmodels/PROVENANCE.md` for the source, the
reference-generation flow, and which upstream NFsim tests aren't
applicable to RM.

Usage:
  python3 harness/basicmodels.py                   # all 29 models
  python3 harness/basicmodels.py --batch 1..20     # contiguous slice
  python3 harness/basicmodels.py --reps 10         # more RM reps
  python3 harness/basicmodels.py r07 r25 r31       # specific models

  # Any flag accepted by benchmark_full.py works (the wrapper passes
  # them through). Useful: --output PATH, --reps N.

The default report path is build/basicmodels_report.md (untracked).
"""

from __future__ import annotations

import argparse
import os
import sys

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
REF_DIR = os.path.join(REPO_ROOT, "tests", "reference", "basicmodels")
DEFAULT_REPORT = os.path.join(REPO_ROOT, "build", "basicmodels_report.md")


def parse_batch(spec: str) -> list[str]:
    """`1..20` → ['r01', ..., 'r20']; `5` → ['r05']."""
    if ".." in spec:
        lo, hi = spec.split("..", 1)
        return [f"r{n:02d}" for n in range(int(lo), int(hi) + 1)]
    return [f"r{int(spec):02d}"]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("models", nargs="*", help="Specific ids like r05 r12; omit for all")
    p.add_argument("--batch", help="Range like '1..20'")
    p.add_argument("--reps", type=int, help="RM replicate count (default 10)")
    p.add_argument(
        "--output",
        default=DEFAULT_REPORT,
        help=f"Report path (default {DEFAULT_REPORT})",
    )
    args, passthrough = p.parse_known_args()

    # Resolve model list. CLI positional wins; otherwise --batch; else all.
    if args.models:
        models = args.models
    elif args.batch:
        models = parse_batch(args.batch)
    else:
        models = []  # benchmark_full.main() reads sim_params.tsv when no models given

    # Build benchmark_full's argv: --ref-dir + (models) + passthrough
    bf_args = ["--ref-dir", REF_DIR, "--output", args.output]
    if args.reps is not None:
        bf_args.extend(["--reps", str(args.reps)])
    bf_args.extend(passthrough)
    bf_args.extend(models)

    # Hand off to benchmark_full.main() with constructed argv. We can't
    # exec it as a subprocess cleanly because benchmark_full reads sys.argv
    # directly; mutate sys.argv and call main() in-process.
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import benchmark_full

    sys.argv = ["benchmark_full.py"] + bf_args
    benchmark_full.main()


if __name__ == "__main__":
    main()
