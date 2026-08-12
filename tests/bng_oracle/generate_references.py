#!/usr/bin/env python3
"""Regenerate the BNG2 references and XML for the tests/bng_oracle suite.

This suite exists because its models are ones where NFsim is *wrong* — they
pin behaviour NFsim does not have — so the usual NFsim ensemble cannot be the
reference.  BioNetGen's own network expansion is, exactly as it is for the
tables in the issues these models come from.

Two reference kinds, chosen per model and recorded in MODELS below:

  ode  `generate_network()` + `simulate(method=>"ode")`.  The default.  Valid
       when the model's molecule counts are large enough that the mass-action
       ODE is the mean of the SSA — which for a catalytic rule with a conserved
       catalyst is exact, since the substrate decay is pseudo-first-order.

  ssa  `generate_network()` + an ensemble of `simulate(method=>"ssa")` reps.
       Needed when the model deliberately runs at a handful of complexes, where
       the ODE's `[X]^2` differs from the SSA's `N(N-1)` and the ODE is
       therefore NOT the mean.  `context_sampler` is that case by construction
       (see its header) and would read 4/3 fast against an ODE reference for
       reasons that have nothing to do with what it tests.  `context_nary` is
       the other case: its rate is a product of two pools that deplete together
       and a catalyst population that fluctuates, so the ODE of the means is
       not the mean either.

Usage:
    BNG2=~/path/to/BNG2.pl python3 tests/bng_oracle/generate_references.py
    ... --model context_sampler      # just one
"""

from __future__ import annotations

import argparse
import concurrent.futures
import glob
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile

SUITE_DIR = os.path.dirname(os.path.abspath(__file__))
MODELS_DIR = os.path.join(SUITE_DIR, "models")
XML_DIR = os.path.join(SUITE_DIR, "xml")
REF_DIR = os.path.join(SUITE_DIR, "reference")

BNG2 = os.environ.get("BNG2", "")

# t_end / n_steps here are the same values check.py drives RuleMonkey with.
MODELS = {
    "context_symmetry": {"kind": "ode", "t_end": 2000, "n_steps": 2},
    "context_sampler": {"kind": "ssa", "t_end": 500, "n_steps": 5, "reps": 2000},
    "context_nary": {"kind": "ssa", "t_end": 100, "n_steps": 2, "reps": 600},
}


def actions_block(body: str) -> str:
    return "begin actions\n" + body + "end actions\n"


def with_actions(text: str, body: str) -> str:
    return re.sub(
        r"begin\s+actions.*?end\s+actions\s*",
        actions_block(body),
        text,
        flags=re.DOTALL,
    )


def run_bng2(bngl_text: str, model: str, timeout: int = 300):
    """Run BNG2 on `bngl_text` in a scratch dir; return (tmpdir contents dir)."""
    tmpdir = tempfile.mkdtemp()
    path = os.path.join(tmpdir, f"{model}.bngl")
    with open(path, "w") as f:
        f.write(bngl_text)
    proc = subprocess.run(
        ["perl", BNG2, f"{model}.bngl"],
        cwd=tmpdir,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout[-4000:] + proc.stderr[-4000:])
        shutil.rmtree(tmpdir, ignore_errors=True)
        raise SystemExit(f"BNG2 failed on {model}")
    return tmpdir


def parse_gdat(path: str):
    header, rows = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                if header is None:
                    header = line.lstrip("#").split()
                continue
            rows.append([float(v) for v in line.split()])
    return header, rows


def write_tsv(path: str, header, rows):
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(f"{v:.10g}" for v in r) + "\n")


def emit_xml(model: str, text: str):
    tmpdir = run_bng2(with_actions(text, "writeXML({overwrite=>1})\n"), model)
    src = os.path.join(tmpdir, f"{model}.xml")
    os.makedirs(XML_DIR, exist_ok=True)
    shutil.copyfile(src, os.path.join(XML_DIR, f"{model}.xml"))
    shutil.rmtree(tmpdir, ignore_errors=True)
    print(f"  xml/{model}.xml")


def emit_ode(model: str, text: str, spec: dict):
    body = (
        "generate_network({overwrite=>1})\n"
        f'simulate({{method=>"ode",t_end=>{spec["t_end"]},n_steps=>{spec["n_steps"]}}})\n'
    )
    tmpdir = run_bng2(with_actions(text, body), model)
    header, rows = parse_gdat(os.path.join(tmpdir, f"{model}.gdat"))
    shutil.rmtree(tmpdir, ignore_errors=True)
    os.makedirs(REF_DIR, exist_ok=True)
    write_tsv(os.path.join(REF_DIR, f"{model}.mean.tsv"), header, rows)
    # A deterministic reference carries no sampling error of its own, so its
    # SEM column is genuinely zero — check.py then scores purely against
    # RuleMonkey's own ensemble spread.
    write_tsv(
        os.path.join(REF_DIR, f"{model}.sem.tsv"),
        header,
        [[r[0]] + [0.0] * (len(r) - 1) for r in rows],
    )
    print(f"  reference/{model}.mean.tsv  (ode, {len(rows)} rows)")


def _ssa_rep(args):
    model, text, t_end, n_steps, seed = args
    body = (
        "generate_network({overwrite=>1})\n"
        f'simulate({{method=>"ssa",t_end=>{t_end},n_steps=>{n_steps},seed=>{seed}}})\n'
    )
    try:
        tmpdir = run_bng2(with_actions(text, body), model, timeout=120)
    except Exception:
        return None
    try:
        gdats = glob.glob(os.path.join(tmpdir, "*.gdat"))
        if not gdats:
            return None
        return parse_gdat(gdats[0])
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def emit_ssa(model: str, text: str, spec: dict):
    reps = spec["reps"]
    jobs = [(model, text, spec["t_end"], spec["n_steps"], seed) for seed in range(1, reps + 1)]
    header, runs = None, []
    n_workers = max(1, (os.cpu_count() or 4) - 1)
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as ex:
        for i, res in enumerate(ex.map(_ssa_rep, jobs, chunksize=4)):
            if res is None:
                continue
            h, rows = res
            header = header or h
            runs.append(rows)
            if (i + 1) % 200 == 0:
                print(f"    {i + 1}/{reps} reps", flush=True)
    if not runs:
        raise SystemExit(f"no SSA reps completed for {model}")

    n_rows = min(len(r) for r in runs)
    n_cols = len(header)
    mean_rows, std_rows = [], []
    for ti in range(n_rows):
        mrow = [runs[0][ti][0]]
        srow = [runs[0][ti][0]]
        for ci in range(1, n_cols):
            vals = [r[ti][ci] for r in runs]
            mrow.append(statistics.mean(vals))
            # SEM, not the per-rep spread: check.py combines it in quadrature
            # with RuleMonkey's own SEM, so the file has to carry the error on
            # the MEAN or the reference's own noise would be double counted.
            srow.append(statistics.stdev(vals) / math.sqrt(len(vals)) if len(vals) > 1 else 0.0)
        mean_rows.append(mrow)
        std_rows.append(srow)
    os.makedirs(REF_DIR, exist_ok=True)
    write_tsv(os.path.join(REF_DIR, f"{model}.mean.tsv"), header, mean_rows)
    write_tsv(os.path.join(REF_DIR, f"{model}.sem.tsv"), header, std_rows)
    print(f"  reference/{model}.mean.tsv  (ssa, {len(runs)} reps)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", action="append", help="regenerate only this model")
    args = ap.parse_args()

    if not BNG2 or not os.path.exists(os.path.expanduser(BNG2)):
        sys.stderr.write("set BNG2 to the path of BNG2.pl\n")
        return 2

    wanted = args.model or sorted(MODELS)
    for model in wanted:
        spec = MODELS[model]
        print(f"{model}:")
        with open(os.path.join(MODELS_DIR, f"{model}.bngl")) as f:
            text = f.read()
        emit_xml(model, text)
        if spec["kind"] == "ode":
            emit_ode(model, text, spec)
        else:
            emit_ssa(model, text, spec)
    return 0


if __name__ == "__main__":
    BNG2 = os.path.expanduser(BNG2)
    sys.exit(main())
