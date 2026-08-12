#!/usr/bin/env python3
"""Score RuleMonkey against the committed BNG2 references in this suite.

Runs `rm_driver` over a fixed seed range per model, averages the trajectories,
and compares each observable at each output time against the reference mean.
The score is a z on the difference of two means,

    z = |rm_mean - ref_mean| / sqrt(rm_sem^2 + ref_sem^2)

so a stochastic reference contributes its own sampling error rather than being
treated as exact.  A deterministic (ODE) reference carries an all-zero SEM
column and the z reduces to RuleMonkey's own ensemble spread.

Exits non-zero on the first model that exceeds the threshold, printing the
worst offender per model either way.

Usage:
    python3 tests/bng_oracle/check.py --driver build/release/rm_driver
"""

from __future__ import annotations

import argparse
import math
import os
import subprocess
import sys

SUITE_DIR = os.path.dirname(os.path.abspath(__file__))
XML_DIR = os.path.join(SUITE_DIR, "xml")
REF_DIR = os.path.join(SUITE_DIR, "reference")

# seeds: enough that the SEM resolves the failure mode this model pins.  For
# context_symmetry that is a factor of 2 or 3 on a pool of 4000, which a
# handful of seeds already separates by tens of sigma; for context_sampler it
# is a 34% shift on 8 molecules, which needs a few hundred.
MODELS = {
    "context_symmetry": {"t_end": 2000, "n_steps": 2, "seeds": 24},
    "context_sampler": {"t_end": 500, "n_steps": 5, "seeds": 400},
    "context_nary": {"t_end": 100, "n_steps": 2, "seeds": 48},
}

Z_THRESHOLD = 6.0


def read_tsv(path: str):
    with open(path) as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip()]
    header = lines[0].split("\t")
    rows = [[float(v) for v in ln.split("\t")] for ln in lines[1:]]
    return header, rows


def run_driver(driver: str, xml: str, t_end, n_steps, seed: int):
    proc = subprocess.run(
        [driver, xml, str(t_end), str(n_steps), str(seed)],
        capture_output=True,
        text=True,
        timeout=300,
    )
    if proc.returncode != 0:
        raise SystemExit(f"rm_driver failed on {xml} seed {seed}:\n{proc.stderr[-2000:]}")
    header, rows = None, []
    for line in proc.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("#"):
            if header is None:
                header = line.lstrip("# ").split("\t")
            continue
        parts = line.split("\t")
        # rm_driver prints a trailing wall-time line with a single field.
        if len(parts) < 2:
            continue
        rows.append([float(v) for v in parts])
    return header, rows


def check_model(driver: str, model: str, spec: dict) -> tuple[bool, str]:
    xml = os.path.join(XML_DIR, f"{model}.xml")
    ref_header, ref_mean = read_tsv(os.path.join(REF_DIR, f"{model}.mean.tsv"))
    _, ref_sem = read_tsv(os.path.join(REF_DIR, f"{model}.sem.tsv"))

    runs = []
    rm_header = None
    for seed in range(1, spec["seeds"] + 1):
        h, rows = run_driver(driver, xml, spec["t_end"], spec["n_steps"], seed)
        rm_header = rm_header or h
        runs.append(rows)

    n_rows = min(len(r) for r in runs)
    # Column i of the reference is named ref_header[i]; find it in RM's output.
    col_of = {name: i for i, name in enumerate(rm_header)}
    worst = (0.0, "")
    for ci, name in enumerate(ref_header):
        if ci == 0 or name not in col_of:
            continue
        rc = col_of[name]
        for ti in range(1, min(n_rows, len(ref_mean))):
            vals = [r[ti][rc] for r in runs]
            n = len(vals)
            mean = sum(vals) / n
            var = sum((v - mean) ** 2 for v in vals) / (n - 1) if n > 1 else 0.0
            rm_sem = math.sqrt(var / n)
            err = math.sqrt(rm_sem**2 + ref_sem[ti][ci] ** 2)
            diff = abs(mean - ref_mean[ti][ci])
            # A pair of exact zeros (an observable neither engine ever moves)
            # has no scale to divide by; call it agreement.
            if err <= 0:
                z = 0.0 if diff <= 1e-9 else float("inf")
            else:
                z = diff / err
            if z > worst[0]:
                worst = (
                    z,
                    f"{name} @ t={ref_mean[ti][0]:g}: rm {mean:.4g} +/- {rm_sem:.3g} "
                    f"vs bng {ref_mean[ti][ci]:.4g} +/- {ref_sem[ti][ci]:.3g}",
                )
    return worst[0] <= Z_THRESHOLD, f"max z {worst[0]:.2f}  [{worst[1]}]"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--driver", required=True)
    ap.add_argument("--model", action="append")
    args = ap.parse_args()

    ok = True
    for model in args.model or sorted(MODELS):
        passed, detail = check_model(args.driver, model, MODELS[model])
        print(f"{'PASS' if passed else 'FAIL'}  {model:20s} {detail}", flush=True)
        ok = ok and passed
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
