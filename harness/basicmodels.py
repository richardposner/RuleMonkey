#!/usr/bin/env python3
"""Exploratory NFsim-parity batch runner for RuleMonkey.

Runs RuleMonkey against the historical NFsim parity test suite — 36
models (`r01.txt`–`r36.txt`, paired with `v01.bngl`–`v36.bngl`)
inherited from `validate/basicModels/` in the nfsim-rm tree. Each
model's BNGL is processed by BNG2.pl to produce an XML model and an
ODE reference trajectory; rm_driver then runs the XML and the
trajectory is compared against the ODE.

This suite has **never been run against RM** before. Treat results as
exploratory — failures may be real RM gaps, BNGL features RM doesn't
support, or harness-wiring issues. Do not gate releases on outcomes.

The harness is intentionally batched: pass `--batch N..M` to limit the
run to a contiguous slice (e.g. `--batch 1..20`). Per-model verdicts
are appended/updated in `tests/models/nfsim_basicmodels/STATUS.tsv` so
incremental runs build up coverage over time.

Usage:
  python3 harness/basicmodels.py --batch 1..20    # first 20 models
  python3 harness/basicmodels.py --batch 21..36   # remaining
  python3 harness/basicmodels.py --batch 5..10    # custom slice
  python3 harness/basicmodels.py --status         # print STATUS.tsv summary

Environment:
  RM_DRIVER  rm_driver binary (default: build/release/rm_driver)
  BNG2       BNG2.pl path     (default: ~/Simulations/BioNetGen-2.9.3/BNG2.pl)
  TIMEOUT    Per-model RM timeout in seconds (default: 60)

Verdicts (column in STATUS.tsv):
  PASS              max relative error to ODE < REL_TOL across observables
  FAIL-divergence   ran cleanly but disagrees with ODE
  FAIL-unsupported  rm_driver exit code 2 (Tier-0 refusal: compartments,
                    Arrhenius, multi-mol Fixed species, ...)
  FAIL-runtime      rm_driver non-zero exit, other reason; see notes
  FAIL-harness      BNG2 didn't produce XML / ODE; not RM's fault
  TIMEOUT           rm_driver exceeded TIMEOUT seconds
"""

from __future__ import annotations

import argparse
import datetime
import glob
import os
import shutil
import subprocess
import sys
import tempfile

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))

RM_DRIVER = os.environ.get(
    "RM_DRIVER",
    os.path.join(REPO_ROOT, "build", "release", "rm_driver"),
)
BNG2 = os.environ.get(
    "BNG2",
    os.path.expanduser("~/Simulations/BioNetGen-2.9.3/BNG2.pl"),
)
TIMEOUT = int(os.environ.get("TIMEOUT", "60"))

SUITE_DIR = os.path.join(REPO_ROOT, "tests", "models", "nfsim_basicmodels")
CACHE_DIR = os.path.join(SUITE_DIR, "_cache")  # gitignored
STATUS_TSV = os.path.join(SUITE_DIR, "STATUS.tsv")

# Comparison tolerance — same 35% threshold used in the legacy NFsim
# parity harness, applied per-observable on max relative error.
REL_TOL = 0.35

STATUS_HEADERS = [
    "model",
    "name",
    "last_run",
    "verdict",
    "rm_s",
    "max_rel_err",
    "worst_obs",
    "notes",
]


# ---------------------------------------------------------------------------
# Configuration / discovery
# ---------------------------------------------------------------------------


def parse_batch(spec: str) -> list[int]:
    """Parse 'N..M' into [N, N+1, ..., M] or 'N' into [N]."""
    if ".." in spec:
        lo, hi = spec.split("..", 1)
        return list(range(int(lo), int(hi) + 1))
    return [int(spec)]


def model_id(num: int) -> str:
    """Zero-padded 2-digit id used in filenames (r01, v01, ...)."""
    return f"{num:02d}"


def load_config(num: int) -> tuple[str, str] | None:
    """Read r{num}.txt; return (model_name, run_options) or None if empty/missing."""
    path = os.path.join(SUITE_DIR, f"r{model_id(num)}.txt")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        lines = [line.strip() for line in f if line.strip()]
    if not lines:
        return None
    name = lines[0]
    opts = lines[1] if len(lines) > 1 else ""
    return name, opts


# ---------------------------------------------------------------------------
# BNG2 caching
# ---------------------------------------------------------------------------


def ensure_bng_outputs(num: int) -> tuple[str, str] | tuple[None, str]:
    """Run BNG2.pl on v{num}.bngl. Outputs cached at _cache/v{num}/.

    Returns (xml_path, ode_gdat_path) on success, or (None, error_msg) on failure.
    """
    mid = model_id(num)
    bngl = os.path.join(SUITE_DIR, f"v{mid}.bngl")
    if not os.path.exists(bngl):
        return None, f"v{mid}.bngl not vendored"

    cache = os.path.join(CACHE_DIR, f"v{mid}")
    xml_path = os.path.join(cache, f"v{mid}.xml")
    ode_path = os.path.join(cache, f"v{mid}_ode.gdat")
    if os.path.exists(xml_path) and os.path.exists(ode_path):
        return xml_path, ode_path

    if not os.path.exists(BNG2):
        return None, f"BNG2 not found at {BNG2}"

    os.makedirs(cache, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        shutil.copy(bngl, tmp)
        try:
            result = subprocess.run(
                ["perl", BNG2, "-outdir", tmp, "-log", os.path.join(tmp, f"v{mid}.bngl")],
                capture_output=True,
                text=True,
                timeout=120,
            )
        except subprocess.TimeoutExpired:
            return None, "BNG2 timed out"
        if result.returncode != 0:
            tail = (result.stderr or result.stdout or "")[-200:].replace("\n", " ")
            return None, f"BNG2 nonzero exit: {tail}"

        xml_files = glob.glob(os.path.join(tmp, "*.xml"))
        if not xml_files:
            return None, "BNG2 produced no XML"
        shutil.copy(xml_files[0], xml_path)

        ode_files = glob.glob(os.path.join(tmp, "*_ode.gdat"))
        if not ode_files:
            return None, "BNG2 produced no ODE output (no simulate_ode action?)"
        shutil.copy(ode_files[0], ode_path)

    return xml_path, ode_path


# ---------------------------------------------------------------------------
# rm_driver invocation + comparison
# ---------------------------------------------------------------------------


def parse_gdat(text: str) -> tuple[list[str], list[list[float]]]:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        return [], []
    headers = lines[0].lstrip("#").strip().split()
    rows: list[list[float]] = []
    for line in lines[1:]:
        if line.startswith("#"):
            continue
        rows.append([float(x) for x in line.split()])
    return headers, rows


def load_gdat_file(path: str) -> tuple[list[str], list[list[float]]]:
    with open(path) as f:
        return parse_gdat(f.read())


def parse_run_options(opts: str) -> tuple[float, int]:
    """Extract t_end and n_steps from r{num}.txt's options line."""
    t_end = 100.0
    n_steps = 100
    toks = opts.split()
    for i, tok in enumerate(toks):
        if tok == "-sim" and i + 1 < len(toks):
            t_end = float(toks[i + 1])
        elif tok == "-oSteps" and i + 1 < len(toks):
            n_steps = int(toks[i + 1])
    return t_end, n_steps


def run_rm(xml_path: str, t_end: float, n_steps: int, seed: int = 42) -> tuple[int, str, str]:
    """Run rm_driver, return (returncode, stdout, stderr)."""
    try:
        result = subprocess.run(
            [RM_DRIVER, xml_path, str(t_end), str(n_steps), str(seed)],
            capture_output=True,
            text=True,
            timeout=TIMEOUT,
        )
        return result.returncode, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return -1, "", "TIMEOUT"


def compare_rm_to_ode(
    rm_h: list[str],
    rm_d: list[list[float]],
    ode_h: list[str],
    ode_d: list[list[float]],
) -> tuple[float, str]:
    """Max relative error across shared observables. Returns (max_err, worst_obs)."""
    if not rm_d or not ode_d:
        return float("inf"), "no-data"

    # Match observable columns by name (skip column 0 = time).
    name_to_ode = {h: i for i, h in enumerate(ode_h)}
    max_err = 0.0
    worst = ""

    # Use the final time-step values for the comparison — sufficient for
    # exploratory pass/fail signal on these small models.
    rm_final = rm_d[-1]
    ode_final = ode_d[-1]

    for j, h in enumerate(rm_h):
        if j == 0 or h not in name_to_ode:
            continue
        rm_v = rm_final[j]
        ode_v = ode_final[name_to_ode[h]]
        denom = max(abs(ode_v), 1.0)  # avoid div-by-zero on rare species
        err = abs(rm_v - ode_v) / denom
        if err > max_err:
            max_err = err
            worst = h
    return max_err, worst


# ---------------------------------------------------------------------------
# STATUS.tsv read/write
# ---------------------------------------------------------------------------


def load_status() -> dict[str, dict[str, str]]:
    """STATUS.tsv → {model_id: {col: value}}."""
    if not os.path.exists(STATUS_TSV):
        return {}
    with open(STATUS_TSV) as f:
        lines = [line.rstrip("\n") for line in f if line.strip()]
    if not lines:
        return {}
    headers = lines[0].split("\t")
    out = {}
    for line in lines[1:]:
        cells = line.split("\t")
        row = {h: cells[i] if i < len(cells) else "" for i, h in enumerate(headers)}
        out[row["model"]] = row
    return out


def write_status(status: dict[str, dict[str, str]]) -> None:
    rows = sorted(status.values(), key=lambda r: r["model"])
    with open(STATUS_TSV, "w") as f:
        f.write("\t".join(STATUS_HEADERS) + "\n")
        for row in rows:
            f.write("\t".join(row.get(h, "") for h in STATUS_HEADERS) + "\n")


# ---------------------------------------------------------------------------
# Per-model run + classification
# ---------------------------------------------------------------------------


def classify_and_run(num: int) -> dict[str, str]:
    """Run one model end-to-end, return a STATUS row."""
    mid = model_id(num)
    cfg = load_config(num)
    name = cfg[0] if cfg else "(no config)"
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    row: dict[str, str] = {
        "model": mid,
        "name": name,
        "last_run": now,
        "verdict": "",
        "rm_s": "",
        "max_rel_err": "",
        "worst_obs": "",
        "notes": "",
    }

    if cfg is None:
        row["verdict"] = "SKIP"
        row["notes"] = "empty config / no r{N}.txt"
        return row

    bng_result = ensure_bng_outputs(num)
    if bng_result[0] is None:
        row["verdict"] = "FAIL-harness"
        row["notes"] = bng_result[1]
        return row
    xml_path, ode_path = bng_result

    t_end, n_steps = parse_run_options(cfg[1])

    import time

    t0 = time.monotonic()
    rc, stdout, stderr = run_rm(xml_path, t_end, n_steps)
    rm_s = time.monotonic() - t0
    row["rm_s"] = f"{rm_s:.2f}"

    if rc == -1:
        row["verdict"] = "TIMEOUT"
        row["notes"] = f"timeout after {TIMEOUT}s"
        return row
    if rc == 2:
        row["verdict"] = "FAIL-unsupported"
        # Pull the most informative line from stderr.
        msg = (stderr or "").splitlines()
        first = next((line for line in msg if line.strip()), "")
        row["notes"] = first[:120]
        return row
    if rc != 0:
        row["verdict"] = "FAIL-runtime"
        msg = (stderr or "").splitlines()
        first = next((line for line in msg if line.strip()), "")
        row["notes"] = first[:120]
        return row

    # Compare to ODE
    rm_h, rm_d = parse_gdat(stdout)
    ode_h, ode_d = load_gdat_file(ode_path)
    max_err, worst = compare_rm_to_ode(rm_h, rm_d, ode_h, ode_d)
    row["max_rel_err"] = f"{max_err:.3f}"
    row["worst_obs"] = worst

    if max_err < REL_TOL:
        row["verdict"] = "PASS"
    else:
        row["verdict"] = "FAIL-divergence"
        row["notes"] = f"final-step rel err {max_err:.2f} > {REL_TOL}"

    return row


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def cmd_status() -> None:
    status = load_status()
    if not status:
        print("STATUS.tsv is empty or missing — run a batch first.")
        return
    by_verdict: dict[str, list[str]] = {}
    for mid, row in sorted(status.items()):
        by_verdict.setdefault(row["verdict"], []).append(mid)
    print(f"{'verdict':<20}  count  models")
    print("-" * 70)
    for verdict in sorted(by_verdict):
        mids = by_verdict[verdict]
        print(f"{verdict:<20}  {len(mids):>5}  {','.join(mids)}")
    print(f"\nTotal: {len(status)} models tracked. See {STATUS_TSV}")


def cmd_batch(nums: list[int]) -> None:
    if not os.path.exists(RM_DRIVER):
        sys.exit(f"rm_driver not built at {RM_DRIVER}; cmake --build --preset release first")
    status = load_status()
    print(f"Running batch: {','.join(model_id(n) for n in nums)}\n")
    print(f"{'model':<6} {'verdict':<18} {'rm_s':>6}  {'rel_err':>8}  notes")
    print("-" * 100)
    for num in nums:
        row = classify_and_run(num)
        status[row["model"]] = row
        print(
            f"r{row['model']:<5} {row['verdict']:<18} {row['rm_s']:>6}  "
            f"{row['max_rel_err']:>8}  {row['notes'][:50]}"
        )
        write_status(status)  # write incrementally so a crash preserves progress

    # Summary
    print("\n" + "-" * 100)
    by_verdict: dict[str, int] = {}
    for num in nums:
        v = status[model_id(num)]["verdict"]
        by_verdict[v] = by_verdict.get(v, 0) + 1
    summary = "  ".join(f"{v}: {n}" for v, n in sorted(by_verdict.items()))
    print(f"Batch summary: {summary}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--batch", help="Range like '1..20' or single id like '5'")
    p.add_argument("--status", action="store_true", help="Print STATUS.tsv summary")
    args = p.parse_args()

    if args.status:
        cmd_status()
        return
    if args.batch:
        cmd_batch(parse_batch(args.batch))
        return
    p.print_help()


if __name__ == "__main__":
    main()
