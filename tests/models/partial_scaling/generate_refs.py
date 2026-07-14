#!/usr/bin/env python3
"""Generate BNG2 reference data for the partial-scaling test models.

For each model in tests/models/partial_scaling/*.bngl this produces, under
tests/models/partial_scaling/reference/:

  net/<model>.net            BNG2 reaction network (generate_network).
  ode/<model>.gdat           BNG2 ODE trajectory (deterministic mean field).
  ssa/<model>.mean.tsv       Ensemble mean over N SSA replicates (seeds 1..N).
  ssa/<model>.std.tsv        Ensemble population std over the same replicates.

BNG2's SSA is a *network-based* Gillespie simulator, so it is an oracle that
is genuinely independent of both RuleMonkey and NFsim (which share the
network-free RuleMonkey-method lineage).  This is the reference the scaled
(Nc > 0) partial-scaling path is validated against, together with the `nfa`
prototype and — where they exist in closed form — analytic moments.  NFsim
is deliberately NOT used here: it has no partial scaling and is not an
oracle for the scaled path.

Usage:
  python3 tests/models/partial_scaling/generate_refs.py            # both models
  python3 tests/models/partial_scaling/generate_refs.py lin2019_toy
  python3 tests/models/partial_scaling/generate_refs.py --reps 500 --workers 8

Environment:
  BNG2   BNG2.pl path (default: ~/Simulations/BioNetGen-2.9.3/BNG2.pl)
"""

from __future__ import annotations

import argparse
import glob
import math
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE = os.path.abspath(os.path.dirname(__file__))
REF_DIR = os.path.join(HERE, "reference")
NET_DIR = os.path.join(REF_DIR, "net")
ODE_DIR = os.path.join(REF_DIR, "ode")
SSA_DIR = os.path.join(REF_DIR, "ssa")

BNG2 = os.environ.get("BNG2", os.path.expanduser("~/Simulations/BioNetGen-2.9.3/BNG2.pl"))
DEFAULT_REPS = 200


def model_body(src_path: str) -> str:
    """Return the model text with the trailing `begin actions` block stripped."""
    text = open(src_path, encoding="utf-8").read()
    idx = text.find("begin actions")
    return text[:idx] if idx >= 0 else text


def run_bng2(bngl_text: str, workdir: str, stem: str) -> None:
    """Write `bngl_text` to workdir/<stem>.bngl and run BNG2 on it."""
    path = os.path.join(workdir, stem + ".bngl")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(bngl_text)
    res = subprocess.run(
        ["perl", BNG2, path],
        cwd=workdir,
        capture_output=True,
        text=True,
        timeout=600,
    )
    if res.returncode != 0:
        raise RuntimeError(f"BNG2 failed for {stem}:\n{res.stdout}\n{res.stderr}")


def read_gdat(path: str) -> tuple[list[str], list[list[float]]]:
    """Parse a BNG .gdat: return (column headers, rows-of-floats)."""
    headers: list[str] = []
    rows: list[list[float]] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.lstrip().startswith("#"):
                headers = line.lstrip("#").split()
                continue
            rows.append([float(x) for x in line.split()])
    return headers, rows


def write_tsv(path: str, headers: list[str], rows: list[list[float]]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(headers) + "\n")
        for row in rows:
            fh.write("\t".join(_fmt(x) for x in row) + "\n")


def _fmt(x: float) -> str:
    if x == int(x):
        return str(int(x))
    return repr(round(x, 6))


def gen_net_and_ode(name: str, src_path: str, t_end: float, n_steps: int) -> None:
    body = model_body(src_path)
    with tempfile.TemporaryDirectory() as wd:
        actions = (
            "begin actions\n"
            "    generate_network({overwrite=>1})\n"
            f'    simulate({{method=>"ode",t_end=>{t_end},n_steps=>{n_steps}}})\n'
            "end actions\n"
        )
        run_bng2(body + actions, wd, name)
        shutil.copyfile(os.path.join(wd, name + ".net"), os.path.join(NET_DIR, name + ".net"))
        # BNG writes the ODE trajectory to <name>.gdat.
        headers, rows = read_gdat(os.path.join(wd, name + ".gdat"))
        write_tsv(os.path.join(ODE_DIR, name + ".gdat"), headers, rows)
    print(f"  [{name}] net + ode written")


def _ssa_rep(args: tuple[str, str, int, float, int]) -> list[list[float]]:
    name, body, seed, t_end, n_steps = args
    with tempfile.TemporaryDirectory() as wd:
        actions = (
            "begin actions\n"
            "    generate_network({overwrite=>1})\n"
            f'    simulate({{method=>"ssa",t_end=>{t_end},n_steps=>{n_steps},seed=>{seed}}})\n'
            "end actions\n"
        )
        run_bng2(body + actions, wd, name)
        _, rows = read_gdat(os.path.join(wd, name + ".gdat"))
    return rows


def gen_ssa_ensemble(
    name: str, src_path: str, t_end: float, n_steps: int, reps: int, workers: int
) -> None:
    body = model_body(src_path)
    tasks = [(name, body, seed, t_end, n_steps) for seed in range(1, reps + 1)]
    reps_rows: list[list[list[float]]] = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = [ex.submit(_ssa_rep, t) for t in tasks]
        for i, fut in enumerate(as_completed(futures), 1):
            reps_rows.append(fut.result())
            if i % 25 == 0 or i == reps:
                print(f"  [{name}] ssa {i}/{reps}")

    # Column headers: reuse the ODE reference already written for this model
    # (same observable order); normalize the leading "#time" to "time".
    headers = _read_tsv_header(os.path.join(ODE_DIR, name + ".gdat"))
    headers[0] = "time"

    n_t = min(len(r) for r in reps_rows)
    n_c = len(headers)
    means: list[list[float]] = []
    stds: list[list[float]] = []
    for ti in range(n_t):
        row_mean = [reps_rows[0][ti][0]]  # time from first rep (grid is shared)
        row_std = [reps_rows[0][ti][0]]
        for ci in range(1, n_c):
            vals = [reps_rows[r][ti][ci] for r in range(len(reps_rows))]
            mu = sum(vals) / len(vals)
            var = sum((v - mu) ** 2 for v in vals) / len(vals)  # population variance
            row_mean.append(mu)
            row_std.append(math.sqrt(var))
        means.append(row_mean)
        stds.append(row_std)

    write_tsv(os.path.join(SSA_DIR, name + ".mean.tsv"), headers, means)
    write_tsv(os.path.join(SSA_DIR, name + ".std.tsv"), headers, stds)
    print(f"  [{name}] ssa mean/std written ({reps} reps)")


def _read_tsv_header(path: str) -> list[str]:
    """First (header) line of a written TSV, split on tabs."""
    with open(path, encoding="utf-8") as fh:
        return fh.readline().rstrip("\n").split("\t")


# (t_end, n_steps) per model — mirror each model's own simulate() action.
MODEL_SPECS = {
    "lin2019_toy": (5.0, 100),
    "birth_death": (10.0, 100),
    "binding": (10.0, 100),
    "homodimer": (10.0, 100),
    "intra_ring": (10.0, 100),  # Phase 3: uni ring close/open, analytic 500/500
    "ab_ring": (10.0, 100),  # Phase 3: disjoint ring close/open + assoc/dissoc, all=200
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("models", nargs="*", help="model stems (default: all)")
    ap.add_argument("--reps", type=int, default=DEFAULT_REPS)
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--skip-ssa", action="store_true", help="net + ode only")
    args = ap.parse_args()

    if not os.path.exists(BNG2):
        print(f"BNG2 not found at {BNG2}; set the BNG2 env var.", file=sys.stderr)
        return 2

    for d in (NET_DIR, ODE_DIR, SSA_DIR):
        os.makedirs(d, exist_ok=True)

    stems = args.models or sorted(
        os.path.splitext(os.path.basename(p))[0] for p in glob.glob(os.path.join(HERE, "*.bngl"))
    )
    for name in stems:
        src = os.path.join(HERE, name + ".bngl")
        if not os.path.exists(src):
            print(f"skip {name}: no {src}", file=sys.stderr)
            continue
        t_end, n_steps = MODEL_SPECS.get(name, (10.0, 100))
        print(f"[{name}] t_end={t_end} n_steps={n_steps}")
        gen_net_and_ode(name, src, t_end, n_steps)
        if not args.skip_ssa:
            gen_ssa_ensemble(name, src, t_end, n_steps, args.reps, args.workers)
    print("done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
