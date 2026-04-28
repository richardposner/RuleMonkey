#!/usr/bin/env python3
"""Generate 100-rep NFsim reference ensembles for the basicmodels suite.

Produces tests/reference/basicmodels/{xml,ensemble,sim_params.tsv,PROVENANCE.md}
analogous to tests/reference/nfsim/ for the corpus suite. Output is the
gold standard against which RuleMonkey is later validated.

Per-model flow:
  1. Run BNG2.pl on tests/models/nfsim_basicmodels/v{NN}.bngl, take its XML.
  2. Translate the run options from r{NN}.txt (NFsim → RM-compatible flags).
  3. Run NFsim 100 reps with translated flags, seeds 1..100.
  4. Aggregate to mean.tsv + std.tsv per model.
  5. Append a row to sim_params.tsv recording what was actually run.

The 31 BNGL files under tests/models/nfsim_basicmodels/ are derived
from NFsim's `validate/basicModels/` regression suite. See
`tests/reference/basicmodels/PROVENANCE.md` for the source, the
upstream NFsim tests not carried over, and the rationale.

Flag translation:
  Pass through to NFsim:  -sim, -oSteps, -cb, -bscb, numeric -gml
  Pass through to RM:     -bscb (default), set_molecule_limit(N) for -gml N
  Strip from both:        -utl 3 (let auto-compute on both sides)
  Strip output-only:      -ss, -rxnlog, -logbuffer, -connect,
                          -trackconnected, -printconnected,
                          -printmoltypes, -printrxncounts

Usage:
  python3 harness/generate_basicmodels_refs.py                  # all
  python3 harness/generate_basicmodels_refs.py r01 r05 r10      # specific ids
  python3 harness/generate_basicmodels_refs.py --workers 4      # parallelism
  python3 harness/generate_basicmodels_refs.py --reps 20        # quick smoke
  python3 harness/generate_basicmodels_refs.py --tint-only      # post-process

Environment:
  NFSIM_BIN  NFsim binary (default: ~/Code/nfsim-rm/build/NFsim — has UTL+1 fix)
  BNG2       BNG2.pl path (default: ~/Simulations/BioNetGen-2.9.3/BNG2.pl)
"""

from __future__ import annotations

import argparse
import datetime
import glob
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SUITE_DIR = os.path.join(REPO_ROOT, "tests", "models", "nfsim_basicmodels")
REF_DIR = os.path.join(REPO_ROOT, "tests", "reference", "basicmodels")
XML_DIR = os.path.join(REF_DIR, "xml")
ENS_DIR = os.path.join(REF_DIR, "ensemble")
REP_DIR = os.path.join(REF_DIR, "replicates")  # gitignored
PARAMS_TSV = os.path.join(REF_DIR, "sim_params.tsv")
PROVENANCE = os.path.join(REF_DIR, "PROVENANCE.md")

NFSIM = os.environ.get("NFSIM_BIN", os.path.expanduser("~/Code/nfsim-rm/build/NFsim"))
BNG2 = os.environ.get("BNG2", os.path.expanduser("~/Simulations/BioNetGen-2.9.3/BNG2.pl"))
N_REPS = 100

# Output-only or -utl override flags to strip from NFsim invocation.
STRIP_FLAGS = {
    "-ss",  # species file output (takes a path argument)
    "-rxnlog",  # reaction firing log JSON (takes a path argument)
    "-logbuffer",  # rxnlog buffer size (takes int)
    "-connect",  # infer connectivity (no argument)
    "-trackconnected",  # rxnlog companion (no argument)
    "-printconnected",  # rxnlog companion (no argument)
    "-printmoltypes",  # output moltypes (no argument)
    "-printrxncounts",  # output rxn counts (no argument)
    "-utl",  # explicit UTL override — let auto-compute on both sides
}
# Flags above that consume one positional arg.
STRIP_FLAGS_WITH_ARG = {"-ss", "-rxnlog", "-logbuffer", "-utl"}


# ---------------------------------------------------------------------------
# Configuration helpers
# ---------------------------------------------------------------------------


def model_id(num: int) -> str:
    return f"r{num:02d}"


def all_model_ids() -> list[str]:
    """Discover testable model ids from the filesystem (r{NN}.txt files)."""
    rids = []
    for path in sorted(glob.glob(os.path.join(SUITE_DIR, "r??.txt"))):
        rids.append(os.path.basename(path)[:-4])
    return rids


def load_config(rid: str) -> tuple[str, str] | None:
    """Read tests/models/nfsim_basicmodels/{rid}.txt → (name, run_options)."""
    path = os.path.join(SUITE_DIR, f"{rid}.txt")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        lines = [line.strip() for line in f if line.strip()]
    if not lines:
        return None
    return lines[0], lines[1] if len(lines) > 1 else ""


def translate_flags(opts: str) -> list[str]:
    """Strip output-only flags + -utl override; return tokens to pass to NFsim."""
    toks = opts.split()
    out: list[str] = []
    i = 0
    while i < len(toks):
        t = toks[i]
        if t in STRIP_FLAGS:
            i += 2 if t in STRIP_FLAGS_WITH_ARG else 1
            continue
        out.append(t)
        i += 1
    # NFsim defaults to -cb off; ensure -bscb (which implies -cb) is present so
    # the gold-standard reference exercises strict BNGL semantics on both sides
    # of every reversible rule, matching RM's default behavior.
    if "-bscb" not in out and "-cb" not in out:
        out.append("-bscb")
    return out


# ---------------------------------------------------------------------------
# BNG2 XML generation
# ---------------------------------------------------------------------------


def generate_xml(rid: str) -> str | None:
    """Run BNG2 -xml on v{NN}.bngl; copy result to XML_DIR/{rid}.xml."""
    num = int(rid[1:])
    bngl = os.path.join(SUITE_DIR, f"v{num:02d}.bngl")
    if not os.path.exists(bngl):
        return None
    out_path = os.path.join(XML_DIR, f"{rid}.xml")
    if os.path.exists(out_path):
        return out_path
    os.makedirs(XML_DIR, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        shutil.copy(bngl, tmp)
        result = subprocess.run(
            ["perl", BNG2, "-outdir", tmp, "-log", "-xml", os.path.basename(bngl)],
            capture_output=True,
            text=True,
            timeout=120,
            cwd=tmp,
        )
        if result.returncode != 0:
            return None
        xml_files = glob.glob(os.path.join(tmp, "*.xml"))
        if not xml_files:
            return None
        shutil.copy(xml_files[0], out_path)
    return out_path


# ---------------------------------------------------------------------------
# NFsim invocation + ensemble aggregation
# ---------------------------------------------------------------------------


def run_nfsim_rep(xml: str, flags: list[str], seed: int, timeout: int = 600) -> str | None:
    """Run one NFsim rep in a temp dir; return .gdat text or None.

    Accepts reps where NFsim produced a valid output file even if the process
    exited nonzero. r31 specifically segfaults on cleanup (rc=-11) but writes
    a complete trajectory first; rejecting that would discard a valid
    reference. Validity check: file exists, has a header line, and at least
    one data row.
    """
    with tempfile.TemporaryDirectory() as tmp:
        out_gdat = os.path.join(tmp, "out.gdat")
        cmd = [NFSIM, "-xml", xml, "-o", out_gdat, "-seed", str(seed)] + flags
        try:
            subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            return None
        if not os.path.exists(out_gdat):
            return None
        with open(out_gdat) as f:
            text = f.read()
        # Validate the trajectory has a header + ≥1 row.
        lines = [line for line in text.splitlines() if line.strip()]
        if len(lines) < 2:
            return None
        return text


def parse_gdat(text: str) -> tuple[list[str], list[list[float]]]:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        return [], []
    headers = lines[0].lstrip("#").strip().split()
    rows: list[list[float]] = []
    for line in lines[1:]:
        if line.startswith("#"):
            continue
        rows.append([float(v) for v in line.split()])
    return headers, rows


def aggregate(reps: list[tuple[list[str], list[list[float]]]]):
    headers = reps[0][0]
    min_times = min(len(r[1]) for r in reps)
    n_obs = len(headers)
    means: list[list[float]] = []
    stds: list[list[float]] = []
    for ti in range(min_times):
        mrow = []
        srow = []
        for ci in range(n_obs):
            vals = [reps[k][1][ti][ci] for k in range(len(reps))]
            mrow.append(statistics.mean(vals))
            srow.append(statistics.stdev(vals) if len(vals) > 1 else 0.0)
        means.append(mrow)
        stds.append(srow)
    return headers, means, stds


def write_tsv(path: str, headers: list[str], rows: list[list[float]]) -> None:
    with open(path, "w") as f:
        f.write("\t".join(headers) + "\n")
        for r in rows:
            f.write("\t".join(f"{v:g}" for v in r) + "\n")


def trapezoidal_integral(times: list[float], values: list[float]) -> float:
    """Trapezoidal rule ∫y dt over the given time grid."""
    n = len(times)
    if n < 2:
        return 0.0
    total = 0.0
    for i in range(n - 1):
        total += 0.5 * (values[i] + values[i + 1]) * (times[i + 1] - times[i])
    return total


def write_tint(rid: str, headers: list[str], reps: list[list[list[float]]]) -> None:
    """Write per-observable time-integral stats: I_mean, I_std, n_reps.

    Format matches what benchmark_full.py's load_tint() expects:
        obs<TAB>I_mean<TAB>I_std<TAB>n_reps
    """
    if not reps:
        return
    # Each rep is a list of [time, obs1, obs2, ...] rows.
    # Compute one I value per rep per observable via trapezoidal integration.
    n_obs = len(headers)
    if n_obs < 2:
        return  # only `time` column; nothing to integrate

    per_obs_I: dict[str, list[float]] = {h: [] for h in headers[1:]}
    for rep_rows in reps:
        if len(rep_rows) < 2:
            continue
        times = [row[0] for row in rep_rows]
        for col in range(1, n_obs):
            values = [row[col] for row in rep_rows]
            per_obs_I[headers[col]].append(trapezoidal_integral(times, values))

    out_path = os.path.join(ENS_DIR, f"{rid}.tint.tsv")
    with open(out_path, "w") as f:
        f.write("obs\tI_mean\tI_std\tn_reps\n")
        for obs_name in headers[1:]:
            integrals = per_obs_I[obs_name]
            if not integrals:
                continue
            n = len(integrals)
            mean = statistics.mean(integrals)
            std = statistics.stdev(integrals) if n > 1 else 0.0
            f.write(f"{obs_name}\t{mean:g}\t{std:g}\t{n}\n")


# ---------------------------------------------------------------------------
# Per-model driver
# ---------------------------------------------------------------------------


def regenerate(rid: str, n_reps: int) -> dict:
    """Returns a result dict for sim_params.tsv + summary."""
    cfg = load_config(rid)
    if cfg is None:
        return {"rid": rid, "status": "no-config", "name": "", "wall_s": 0}
    name, opts = cfg

    xml = generate_xml(rid)
    if xml is None:
        return {"rid": rid, "status": "bng2-failed", "name": name, "wall_s": 0}

    flags = translate_flags(opts)
    # Track the t_end / n_steps explicitly so the harness can read them.
    t_end = "100"
    n_steps = "100"
    for i, t in enumerate(flags):
        if t == "-sim" and i + 1 < len(flags):
            t_end = flags[i + 1]
        elif t == "-oSteps" and i + 1 < len(flags):
            n_steps = flags[i + 1]

    rep_dir = os.path.join(REP_DIR, rid)
    os.makedirs(rep_dir, exist_ok=True)

    parsed_reps: list[tuple[list[str], list[list[float]]]] = []
    t0 = time.monotonic()
    failed = 0
    for i in range(n_reps):
        seed = i + 1
        text = run_nfsim_rep(xml, flags, seed)
        if text is None:
            failed += 1
            continue
        # Save replicate text
        with open(os.path.join(rep_dir, f"rep_{seed:03d}.gdat"), "w") as f:
            f.write(text)
        parsed_reps.append(parse_gdat(text))
        sys.stdout.write(f"\r  {rid}: {len(parsed_reps)}/{n_reps} reps")
        sys.stdout.flush()
    sys.stdout.write("\r" + " " * 60 + "\r")
    sys.stdout.flush()

    if not parsed_reps:
        return {"rid": rid, "status": "all-reps-failed", "name": name, "wall_s": 0}

    headers, means, stds = aggregate(parsed_reps)
    os.makedirs(ENS_DIR, exist_ok=True)
    write_tsv(os.path.join(ENS_DIR, f"{rid}.mean.tsv"), headers, means)
    write_tsv(os.path.join(ENS_DIR, f"{rid}.std.tsv"), headers, stds)
    # Per-observable trapezoidal-integral stats consumed by benchmark_full's
    # tz_max verdict math.
    write_tint(rid, headers, [r[1] for r in parsed_reps])

    wall = time.monotonic() - t0
    return {
        "rid": rid,
        "status": "ok" if failed == 0 else f"ok ({failed} reps failed)",
        "name": name,
        "wall_s": wall,
        "t_end": t_end,
        "n_steps": n_steps,
        "flags": " ".join(flags),
        "n_reps": len(parsed_reps),
    }


# ---------------------------------------------------------------------------
# Output: sim_params.tsv + PROVENANCE.md
# ---------------------------------------------------------------------------


def write_sim_params(results: list[dict]) -> None:
    rows = []
    for r in results:
        if r["status"] == "excluded":
            continue
        if r["status"].startswith("ok"):
            rows.append(
                "\t".join(
                    [
                        r["rid"],
                        r["name"][:60],
                        r["t_end"],
                        r["n_steps"],
                        r["flags"],
                        f"{r['wall_s']:.1f}",
                        str(r["n_reps"]),
                    ]
                )
            )
    with open(PARAMS_TSV, "w") as f:
        f.write(
            "\t".join(["model", "name", "t_end", "n_steps", "nfsim_flags", "wall_s", "n_reps"])
            + "\n"
        )
        f.write("\n".join(rows) + "\n")


def write_provenance(results: list[dict], n_reps: int) -> None:
    today = datetime.date.today().isoformat()
    body = f"""# basicmodels Reference Data — Provenance

100-replicate NFsim ensemble references for the historical NFsim parity
suite (33 testable models). Three models from the original NFsim
test set were removed because they don't apply to RM:

- **r27, r28** — used the BNGL `population` keyword (hybrid
  particle-population SSA, Hogg 2013). RM has no equivalent and now
  refuses such models at Tier-0 (`MoleculeType@population` error). See
  `cpp/rulemonkey/simulator.cpp:scan_unsupported`.
- **r36** — tested NFsim's `-gml auto` fallback (issue #17). RM only
  honors numeric `set_molecule_limit`; the `auto` token is a
  non-feature, not a bug.

**Last regeneration:** {today}
**NFsim binary:** `{NFSIM}` (has UTL+1 fix)
**BNG2:** `{BNG2}`
**Reps per model:** {n_reps}

## Flag translation policy

The `r{{NN}}.txt` config files specify NFsim run options that include a
mix of trajectory-affecting flags and output-only / debug flags. The
generator translates them so the gold-standard reference reflects what
the trajectory should look like, not which output files NFsim happens
to produce. See `harness/generate_basicmodels_refs.py:STRIP_FLAGS` for
the canonical list.

Notably, `-utl 3` (used by r21–r26 in the original config) is dropped
on both sides — both NFsim (with the UTL+1 patch) and RM auto-compute
`max_pattern_size + 1`, which is the correct value.

## Per-model summary

See `sim_params.tsv` for the actual flags used per model and the wall
time of each NFsim ref-gen run.
"""
    with open(PROVENANCE, "w") as f:
        f.write(body)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _worker(rid: str, n_reps: int) -> dict:
    return regenerate(rid, n_reps)


def regenerate_tint_only(rid: str) -> dict:
    """Recompute tint.tsv from existing on-disk replicates without re-running NFsim."""
    rep_dir = os.path.join(REP_DIR, rid)
    if not os.path.isdir(rep_dir):
        return {"rid": rid, "status": "no-replicates", "name": "", "wall_s": 0}
    rep_files = sorted(glob.glob(os.path.join(rep_dir, "rep_*.gdat")))
    parsed: list[tuple[list[str], list[list[float]]]] = []
    headers: list[str] | None = None
    for rp in rep_files:
        with open(rp) as f:
            text = f.read()
        h, rows = parse_gdat(text)
        if not h:
            continue
        if headers is None:
            headers = h
        elif h != headers:
            continue
        parsed.append((h, rows))
    if not parsed or headers is None:
        return {"rid": rid, "status": "no-valid-reps", "name": "", "wall_s": 0}
    write_tint(rid, headers, [p_[1] for p_ in parsed])
    return {"rid": rid, "status": f"ok ({len(parsed)} reps)", "name": "", "wall_s": 0}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("models", nargs="*", help="Specific ids like r01 r05; omit for all 36")
    p.add_argument("--reps", type=int, default=N_REPS, help="Replicates per model")
    p.add_argument("--workers", type=int, default=4, help="Parallel models")
    p.add_argument(
        "--tint-only",
        action="store_true",
        help="Recompute *.tint.tsv from existing on-disk replicates (no NFsim runs)",
    )
    args = p.parse_args()

    if args.tint_only:
        rids = args.models if args.models else all_model_ids()
        for rid in rids:
            r = regenerate_tint_only(rid)
            print(f"  {rid:<5} {r['status']}")
        return

    if not os.path.exists(NFSIM):
        sys.exit(f"NFSIM binary not found at {NFSIM} (set NFSIM_BIN env var)")
    if not os.path.exists(BNG2):
        sys.exit(f"BNG2 not found at {BNG2} (set BNG2 env var)")

    rids = args.models if args.models else all_model_ids()
    print(f"Generating {len(rids)} models × {args.reps} reps with {args.workers} workers\n")

    os.makedirs(REF_DIR, exist_ok=True)
    results: list[dict] = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futures = {ex.submit(_worker, rid, args.reps): rid for rid in rids}
        for fut in as_completed(futures):
            rid = futures[fut]
            try:
                r = fut.result()
            except Exception as e:
                r = {"rid": rid, "status": f"crashed: {e}", "name": "", "wall_s": 0}
            results.append(r)
            wall = r.get("wall_s", 0)
            status = r["status"]
            print(f"  {rid:<5} {status:<22} {wall:7.1f}s  {r.get('name', '')[:50]}")

    results.sort(key=lambda r: r["rid"])
    write_sim_params(results)
    write_provenance(results, args.reps)

    print(f"\nWrote {PARAMS_TSV}")
    print(f"Wrote {PROVENANCE}")
    print(
        f"\nDone: {sum(1 for r in results if r['status'].startswith('ok'))} ok / "
        f"{sum(1 for r in results if r['status'] == 'excluded')} excluded / "
        f"{sum(1 for r in results if r['status'] not in ('ok',) and not r['status'].startswith('ok') and r['status'] != 'excluded')} failed"
    )


if __name__ == "__main__":
    main()
