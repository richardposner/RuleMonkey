#!/usr/bin/env python3
"""Comprehensive RuleMonkey benchmark: correctness + efficiency.

Runs every model with N replicate seeds, computes:
  - z_mean:    max |RM_mean - NFsim_mean| / (NFsim_std / sqrt(n_reps))
  - std_ratio: max RM_std / NFsim_std  (variance consistency check)
  - Timing:    RM wall time (mean of N reps), phase breakdown
  - Throughput: events/sec

Outputs a markdown report incrementally (results appear as models finish).
Models are run fastest-first so you see progress quickly.

## Tiers (bug-sprint feedback loop)

  --tier smoke  : 8 fast, diverse models, 2 reps each  (~2 min)
                  use while iterating on a fix; catches catastrophic
                  regressions only — z verdicts are noisy at 2 reps
  --tier guard  : 29 regression-guard models, 3 reps   (~15 min)
                  run before every commit — must stay green
  --tier full   : all 71 models, 10 reps               (~3 h)
                  run before declaring a bug fixed / before merge

Explicit --reps overrides the tier default. Positional model args override
the tier model list.  --adaptive inflates reps per model so that noisy
models (high noise floor) get enough reps for T ≈ 5.0, while quiet
models stay at the base rep count.

Usage:
  python3 harness/benchmark_full.py                        # all 71 models, 10 reps
  python3 harness/benchmark_full.py --tier smoke           # 8 models, 2 reps
  python3 harness/benchmark_full.py --tier guard           # 29 models, 3 reps
  python3 harness/benchmark_full.py --tier full --reps 5   # all, 5 reps
  python3 harness/benchmark_full.py --reps 3               # all, 3 reps
  python3 harness/benchmark_full.py --adaptive             # all, 10+ reps per model
  python3 harness/benchmark_full.py tcr ensemble           # specific models
  python3 harness/benchmark_full.py --output report.md     # custom output

Tier membership is defined inline in the SMOKE_MODELS / GUARD_MODELS
lists below — those are the canonical source.
"""

import math
import os
import re
import statistics
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
RM_DRIVER = os.environ.get(
    "RM_DRIVER",
    os.path.join(REPO_ROOT, "build", "release", "rm_driver"),
)
# Reference data root. Default = corpus suite; --ref-dir overrides at runtime.
# Other ref-derived paths are recomputed by set_ref_dir() when --ref-dir is given.
REF_DIR = os.path.join(REPO_ROOT, "tests", "reference", "nfsim")
XML_DIR = os.path.join(REF_DIR, "xml")
PARAMS = os.path.join(REF_DIR, "sim_params.tsv")
ENS_DIR = os.path.join(REF_DIR, "ensemble")
SUMMARY = os.path.join(REF_DIR, "summary.tsv")
NOISE_FLOOR_FILE = os.path.join(REF_DIR, "noise_floor.tsv")
TIMEOUT_FILE = os.path.join(REPO_ROOT, "harness", "model_timeouts.tsv")


def set_ref_dir(new_ref_dir: str) -> None:
    """Re-point all reference-relative paths at a different corpus root.

    Used by basicmodels.py (and any future per-corpus harness) to share
    benchmark_full's verdict logic without code duplication.
    """
    global REF_DIR, XML_DIR, PARAMS, ENS_DIR, SUMMARY, NOISE_FLOOR_FILE
    REF_DIR = os.path.abspath(new_ref_dir)
    XML_DIR = os.path.join(REF_DIR, "xml")
    PARAMS = os.path.join(REF_DIR, "sim_params.tsv")
    ENS_DIR = os.path.join(REF_DIR, "ensemble")
    SUMMARY = os.path.join(REF_DIR, "summary.tsv")
    NOISE_FLOOR_FILE = os.path.join(REF_DIR, "noise_floor.tsv")


DEFAULT_REPS = 10
ADAPTIVE_MAX_REPS = 100  # cap for adaptive per-model reps

# Rep-level parallelism: reps within a single model run concurrently.
# Defaults to min(n_reps, n_cores).  Parallel reps share cache/RAM, so
# per-rep wall inflates modestly (~5-15% vs pure serial on a 6-core
# MBP) but total model wall scales ~linearly with workers.  Override
# via MAX_PARALLEL_REPS env var.
MAX_PARALLEL_REPS = int(os.environ.get("MAX_PARALLEL_REPS", os.cpu_count() or 4))

# Screen (primary): max-z over all (time, observable) pairs. Fast, noisy on
# rare-event Size_N observables. Unchanged from historical behavior; kept as
# an early-warning "something is off" flag.
SCREEN_Z_THRESHOLD = 5.0

# Verdict (secondary): max time-integrated-z over observables, compared to a
# per-model noise floor derived from self-split of the NFsim replicates.
# T_model = max(VERDICT_TZ_FLOOR, VERDICT_TZ_MARGIN * tz_p99) where p99 comes
# from tests/reference/nfsim/noise_floor.tsv.
VERDICT_TZ_FLOOR = 5.0
VERDICT_TZ_MARGIN = 1.2
VERDICT_TZ_DEFAULT = 5.0  # fallback when no calibration row is present

STD_RATIO_WARN = 2.0  # flag if RM_std/NFsim_std outside [1/this, this]
ZERO_EPS = 1e-9
TINT_TOL_REL = 1e-6


# ---------------------------------------------------------------------------
# Tiered model sets — canonical definitions.
# ---------------------------------------------------------------------------

# Smoke tier: fastest diverse coverage. Every major code path touched by at
# least one model. Used for tight-loop iteration while developing a fix.
SMOKE_MODELS = [
    "Tutorial_Example",  # multi-mol observables (UTL)
    "A_plus_A",  # simplest bimolecular
    "BLBR",  # cooperativity + rings
    "isingspin_localfcn",  # local functions
    "e1",  # enzyme kinetics scaling baseline
    "ANx_noActivity",  # ANx without function rate laws (baseline)
    "fceri_ji",  # multi-mol Species + -bscb
    "nfsim_ring_closure_polymer",  # ring closure
]

# Regression guard tier: the models that historically catch regressions.
# Any fix must keep all of these green. Run this tier before every commit.
GUARD_MODELS = [
    # Function rate laws
    "AN",
    "ANx",
    "ANx_noActivity",
    # Enzyme kinetics scaling
    "e1",
    "e2",
    "e3",
    "e4",
    "e5",
    "e6",
    "e7",
    "e8",
    "e9",
    # Multi-mol observables / UTL
    "Tutorial_Example",
    # Multi-mol Species + -bscb
    "fceri_ji",
    # Cooperativity & rings
    "BLBR",
    "bench_blbr_cooperativity_posner2004",
    "bench_blbr_cooperativity_posner2004_rings",
    # Local functions
    "isingspin_localfcn",
    # Ring closure
    "nfsim_ring_closure_polymer",
    "A_plus_B_rings",
    "rm_tlbr_rings",
    # Large rule sets
    "tcr",
    "tcr_gen20ind9",
    # Large complexes
    "egfr_net",
    "egfr_nf_iter5p12h10",
    # High event count stress
    "t3",
    "machine",
    "ensemble",
    # Disjoint transphosphorylation (SSA-validated reference)
    "toy_jim",
]

TIER_DEFAULT_REPS = {
    "smoke": 2,
    "guard": 3,
    "full": DEFAULT_REPS,
}

TIER_MODELS = {
    "smoke": SMOKE_MODELS,
    "guard": GUARD_MODELS,
    "full": None,  # None = all models in sim_params
}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def load_sim_params():
    params = {}
    with open(PARAMS) as f:
        for line in f:
            parts = line.strip().split("\t")
            if parts[0].startswith("#") or parts[0] == "model":
                continue
            if len(parts) >= 4:
                flags = parts[5].strip() if len(parts) >= 6 else ""
                # Extract RM-relevant flags from NFsim flags
                rm_flags = []
                if "-bscb" in flags:
                    rm_flags.append("-bscb")
                params[parts[0]] = (parts[2], parts[3], rm_flags)
    return params


def load_nfsim_wall_times():
    """Load NFsim mean wall times from summary.tsv."""
    times = {}
    if not os.path.exists(SUMMARY):
        return times
    with open(SUMMARY) as f:
        for line in f:
            parts = line.strip().split("\t")
            if parts[0].startswith("#") or parts[0] == "model":
                continue
            if len(parts) >= 8:
                times[parts[0]] = float(parts[6])  # wall_mean_s
    return times


def load_timeouts():
    timeouts = {}
    if not os.path.exists(TIMEOUT_FILE):
        return timeouts
    with open(TIMEOUT_FILE) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("model"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                timeouts[parts[0]] = int(parts[1])
    return timeouts


def read_tsv(path):
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip()]
    headers = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        rows.append([float(v) for v in line.split("\t")])
    return headers, rows


# Sample sizes at which `noise_floor.tsv` is calibrated.  Must match
# the calibration sample sizes used to produce the committed tsv
# (see tests/reference/nfsim/PROVENANCE.md).
CALIBRATED_N = (2, 3, 10)


def load_noise_floor():
    """Load per-model noise-floor calibration from noise_floor.tsv.

    The file has one row per model with tz_p95/p99/max and mz_p95/p99/max
    computed at each sample size in CALIBRATED_N. Returns a dict keyed by
    model name whose values are the raw row dicts (string values); use
    `nf_stat(row, key, n_reps)` to read a numeric metric at the sample
    size closest to `n_reps`.
    """
    nf = {}
    if not os.path.exists(NOISE_FLOOR_FILE):
        return nf
    with open(NOISE_FLOOR_FILE) as f:
        header = None
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
                continue
            row = dict(zip(header, parts, strict=False))
            if "model" in row:
                nf[row["model"]] = row
    return nf


def nf_stat(row, metric, n_reps):
    """Return float value of `metric` calibrated at the sample size in
    CALIBRATED_N that best matches `n_reps`. `metric` is e.g. 'tz_p99'.

    Prefers an exact n_reps match; otherwise falls back to the nearest
    calibrated sample size. Returns None if the row is missing the column.
    """
    if row is None:
        return None
    if n_reps in CALIBRATED_N:
        key = f"{metric}_n{n_reps}"
    else:
        nearest = min(CALIBRATED_N, key=lambda x: abs(x - n_reps))
        key = f"{metric}_n{nearest}"
    if key not in row:
        return None
    try:
        return float(row[key])
    except (TypeError, ValueError):
        return None


def load_tint(model):
    """Load NFsim per-observable time-integral stats for `model`.

    Returns {obs_name: (I_mean, I_std, n_reps)}, or None if missing.
    """
    path = os.path.join(ENS_DIR, f"{model}.tint.tsv")
    if not os.path.exists(path):
        return None
    tint = {}
    with open(path) as f:
        header = None
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
                continue
            try:
                tint[parts[0]] = (float(parts[1]), float(parts[2]), int(parts[3]))
            except (IndexError, ValueError):
                continue
    return tint


def t_model(nf_row, n_reps):
    """Derive the per-model verdict threshold.

    T_model = max(VERDICT_TZ_FLOOR, VERDICT_TZ_MARGIN * tz_p99),
    where tz_p99 is the 99th percentile of `tz_max` under the null, measured
    by self-splitting the 100-rep NFsim reference at a sample size matching
    `n_reps` (via nf_stat). Returns VERDICT_TZ_DEFAULT when no calibration
    row is available.
    """
    if nf_row is None:
        return VERDICT_TZ_DEFAULT
    p99 = nf_stat(nf_row, "tz_p99", n_reps)
    if p99 is None:
        return VERDICT_TZ_DEFAULT
    return max(VERDICT_TZ_FLOOR, VERDICT_TZ_MARGIN * p99)


def adaptive_reps(nf_row, base_reps):
    """Compute minimum reps so that T_model ≈ VERDICT_TZ_FLOOR.

    The noise floor tz_p99 scales as ~1/sqrt(n). Given tz_p99 at n=10,
    find the smallest n >= base_reps such that
    VERDICT_TZ_MARGIN * tz_p99 * sqrt(10/n) <= VERDICT_TZ_FLOOR.
    """
    if nf_row is None:
        return base_reps
    p99_n10 = nf_stat(nf_row, "tz_p99", 10)
    if p99_n10 is None:
        return base_reps
    target = VERDICT_TZ_FLOOR / VERDICT_TZ_MARGIN  # tz_p99 must be <= this
    if p99_n10 <= target:
        return base_reps
    import math

    n = math.ceil(10 * (p99_n10 / target) ** 2)
    return max(base_reps, min(n, ADAPTIVE_MAX_REPS))


def parse_gdat(text):
    lines = [l.strip() for l in text.strip().split("\n") if l.strip()]
    headers = lines[0].lstrip("#").strip().split("\t")
    rows = []
    for line in lines[1:]:
        rows.append([float(v) for v in line.split("\t")])
    return headers, rows


def parse_timing(stderr_text):
    """Extract timing info from RM stderr output."""
    info = {
        "events": 0,
        "null": 0,
        "total_s": 0,
        "select_pct": 0,
        "fire_pct": 0,
        "obs_pct": 0,
        "update_pct": 0,
    }
    m = re.search(r"events=(\d+)", stderr_text)
    if m:
        info["events"] = int(m.group(1))
    m = re.search(r"null=(\d+)", stderr_text)
    if m:
        info["null"] = int(m.group(1))
    m = re.search(r"total=([\d.]+)s", stderr_text)
    if m:
        info["total_s"] = float(m.group(1))
    m = re.search(r"select_reactants:\s*[\d.]+s\s*\(([\d.]+)%\)", stderr_text)
    if m:
        info["select_pct"] = float(m.group(1))
    m = re.search(r"fire_rule:\s*[\d.]+s\s*\(([\d.]+)%\)", stderr_text)
    if m:
        info["fire_pct"] = float(m.group(1))
    m = re.search(r"observables:\s*[\d.]+s\s*\(([\d.]+)%\)", stderr_text)
    if m:
        info["obs_pct"] = float(m.group(1))
    m = re.search(r"incr_update:\s*[\d.]+s\s*\(([\d.]+)%\)", stderr_text)
    if m:
        info["update_pct"] = float(m.group(1))
    return info


# ---------------------------------------------------------------------------
# Run model
# ---------------------------------------------------------------------------


def run_one_rep(xml, t_end, n_steps, seed, timeout=None, rm_flags=None):
    """Run a single RM rep. Returns (headers, rows, stderr_text, wall_s) or None.

    None signals a terminal problem for this rep (timeout, driver crash,
    missing binary, or unparseable output). Anything outside the narrow
    list below is left to propagate so genuine bugs are visible rather
    than silently counted as a timeout.
    """
    t0 = time.monotonic()
    try:
        cmd = [RM_DRIVER, xml, t_end, n_steps, str(seed)]
        if rm_flags:
            cmd.extend(rm_flags)
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        wall_s = time.monotonic() - t0
        if result.returncode != 0 or not result.stdout.strip():
            print(
                f"  [run_one_rep] {os.path.basename(xml)} seed={seed} rc={result.returncode}",
                file=sys.stderr,
            )
            return None
        headers, rows = parse_gdat(result.stdout)
        return (headers, rows, result.stderr, wall_s)
    except subprocess.TimeoutExpired:
        return None
    except (FileNotFoundError, ValueError) as e:
        print(
            f"  [run_one_rep] {os.path.basename(xml)} seed={seed} {type(e).__name__}: {e}",
            file=sys.stderr,
        )
        return None


# ---------------------------------------------------------------------------
# Correctness analysis
# ---------------------------------------------------------------------------


def analyze_correctness(model, all_reps, n_reps, nf_row=None):
    """Compare RM ensemble (n_reps trajectories) against NFsim reference.

    Computes two metrics:
      * screen_z (historical primary): max |RM_mean - NF_mean| /
        (NF_std / sqrt(n_reps)) over all (time, observable) pairs. Fast
        and sensitive, but statistically noisy on models with many
        rare-event Size_N observables.
      * tz_max (verdict): max over observables of the z-score of the
        per-rep trapezoidal time integral, computed against precomputed
        NFsim stats in `ensemble/{model}.tint.tsv`. Compared against a
        per-model noise-floor threshold derived from `nf_row`.

    Returns dict with: screen_z, worst_screen_obs, tz_max, worst_tz_obs,
                       T_model, std_ratio, worst_obs_var, verdict, n_compared.
    """
    mean_path = os.path.join(ENS_DIR, f"{model}.mean.tsv")
    std_path = os.path.join(ENS_DIR, f"{model}.std.tsv")

    if not os.path.exists(mean_path) or not os.path.exists(std_path):
        return {"verdict": "NO_REF"}

    nf_headers, nf_mean = read_tsv(mean_path)
    _, nf_std = read_tsv(std_path)

    rm_headers = all_reps[0][0]

    # Column mapping: RM header -> NF column index
    col_map = {}
    for ri, rh in enumerate(rm_headers):
        for ni, nh in enumerate(nf_headers):
            if rh.strip() == nh.strip():
                col_map[ri] = ni
                break

    # NF time -> row index
    nf_time_idx = {round(row[0], 6): i for i, row in enumerate(nf_mean)}

    # Collect per-timepoint, per-observable values across reps
    # obs_data[obs_name][time_idx] = list of RM values across reps
    obs_data = {}
    time_points = []

    # Use first rep to get time grid
    for rm_row in all_reps[0][1]:
        rm_t = round(rm_row[0], 6)
        if rm_t in nf_time_idx:
            time_points.append(rm_t)

    for rm_ci, nf_ci in col_map.items():
        if rm_ci == 0 or nf_ci == 0:
            continue
        obs_name = rm_headers[rm_ci]
        obs_data[obs_name] = {}
        for tp in time_points:
            obs_data[obs_name][tp] = []

    # Collect values
    for rep_headers, rep_rows, _, _ in all_reps:
        for rm_row in rep_rows:
            rm_t = round(rm_row[0], 6)
            if rm_t not in nf_time_idx:
                continue
            for rm_ci, nf_ci in col_map.items():
                if rm_ci == 0 or nf_ci == 0:
                    continue
                obs_name = rm_headers[rm_ci]
                obs_data[obs_name][rm_t].append(rm_row[rm_ci])

    # Compute z_mean and std_ratio
    max_z_mean = 0
    worst_obs_mean = ""
    max_std_ratio = 0
    worst_obs_var = ""
    n_compared = 0

    for obs_name, tdata in obs_data.items():
        for tp, rm_vals in tdata.items():
            if len(rm_vals) < 1:
                continue
            n_compared += 1
            nf_idx = nf_time_idx[tp]

            nf_ci = None
            for rm_ci, nc in col_map.items():
                if rm_headers[rm_ci] == obs_name:
                    nf_ci = nc
                    break
            if nf_ci is None:
                continue

            nf_m = nf_mean[nf_idx][nf_ci]
            nf_s = nf_std[nf_idx][nf_ci]

            rm_m = statistics.mean(rm_vals)
            rm_s = statistics.stdev(rm_vals) if len(rm_vals) > 1 else 0

            # z-score for means
            if nf_s > 0:
                z = abs(rm_m - nf_m) / (nf_s / math.sqrt(len(rm_vals)))
            elif abs(rm_m - nf_m) < 1e-9:
                z = 0.0
            else:
                z = abs(rm_m - nf_m)

            if z > max_z_mean:
                max_z_mean = z
                worst_obs_mean = obs_name

            # Std ratio (only meaningful when NFsim std > 0 and mean is
            # large enough to be non-trivial)
            if nf_s > 1e-6 and nf_m > 1.0:
                ratio = rm_s / nf_s
                deviation = abs(ratio - 1.0)
                if deviation > max_std_ratio:
                    max_std_ratio = deviation
                    worst_obs_var = f"{obs_name} ({ratio:.2f})"

    if n_compared == 0:
        return {"verdict": "NO_MATCH"}

    # -------- Secondary metric: time-integrated z (tz_max) ------------------
    # Compute each RM rep's trapezoidal integral per observable, then z-score
    # against the NFsim per-observable time-integral stats in the tint file.
    tz_max = 0.0
    worst_tz_obs = ""
    tint = load_tint(model)
    T = t_model(nf_row, n_reps)

    if tint is None:
        tz_max = None
        T = None
        tz_reason = "no-tint"
    else:
        # Per-rep integrals: {obs_name: [I_rep0, I_rep1, ...]}
        # Use the RM time grid (from the first rep) and filter to time points
        # that also appear in the NF time grid, matching the primary metric.
        tp_sorted = sorted(time_points)
        # trapezoidal weights on sorted time grid
        if len(tp_sorted) < 2:
            tz_max = None
            T = None
            tz_reason = "too-few-timepoints"
        else:
            tz_reason = None
            # Index RM rep rows by time for per-rep integration
            per_rep_I = {}  # obs_name -> list of per-rep integrals
            for rm_ci, nf_ci in col_map.items():
                if rm_ci == 0 or nf_ci == 0:
                    continue
                obs_name = rm_headers[rm_ci]
                if obs_name not in tint:
                    continue
                per_rep_I[obs_name] = []

            # Silently-dropped-rep counter: if a rep's time grid misses any
            # shared timepoint, the rep is skipped for that observable. A
            # healthy run should always have len(per_rep_I[obs]) == n_reps;
            # anything less means the verdict for that observable was
            # computed from fewer than the nominal sample count.
            dropped_per_obs = {obs: 0 for obs in per_rep_I}

            for rep_headers, rep_rows, _, _ in all_reps:
                # Build time -> values map for this rep
                t2v = {}
                for rm_row in rep_rows:
                    t2v[round(rm_row[0], 6)] = rm_row

                # Integrate each observable on the shared time grid
                for rm_ci, _nf_ci in col_map.items():
                    if rm_ci == 0:
                        continue
                    obs_name = rm_headers[rm_ci]
                    if obs_name not in per_rep_I:
                        continue
                    # gather values at tp_sorted
                    vals = []
                    ok = True
                    for tp in tp_sorted:
                        if tp not in t2v:
                            ok = False
                            break
                        vals.append(t2v[tp][rm_ci])
                    if not ok:
                        dropped_per_obs[obs_name] += 1
                        continue
                    # trapezoidal integral
                    I = 0.0
                    for i in range(len(tp_sorted) - 1):
                        dt = tp_sorted[i + 1] - tp_sorted[i]
                        I += 0.5 * (vals[i] + vals[i + 1]) * dt
                    per_rep_I[obs_name].append(I)

            # Warn if any observable lost reps.
            any_drops = {o: n for o, n in dropped_per_obs.items() if n > 0}
            if any_drops:
                summary = ", ".join(
                    f"{o}:{n}" for o, n in sorted(any_drops.items(), key=lambda kv: -kv[1])[:5]
                )
                print(
                    f"  [drop] {len(any_drops)} obs lost reps (time-grid mismatch): {summary}",
                    file=sys.stderr,
                )

            for obs_name, I_list in per_rep_I.items():
                if len(I_list) < 1:
                    continue
                nf_mean_I, nf_std_I, n_nf = tint[obs_name]
                rm_mean_I = sum(I_list) / len(I_list)
                if len(I_list) > 1:
                    rm_std_I = statistics.stdev(I_list)
                else:
                    rm_std_I = 0.0
                SE_rm = rm_std_I / math.sqrt(len(I_list))
                SE_nf = nf_std_I / math.sqrt(n_nf) if n_nf > 0 else 0.0
                dm = abs(rm_mean_I - nf_mean_I)
                tol = TINT_TOL_REL * (abs(nf_mean_I) + 1.0)

                # Degenerate cases (both stds zero)
                if nf_std_I < ZERO_EPS and rm_std_I < ZERO_EPS:
                    if dm < tol:
                        continue  # trivially equal
                    # Sustained mismatch between two zero-variance series:
                    # conservation violation or degenerate observable. Always fail.
                    tz_max = float("inf")
                    worst_tz_obs = obs_name + "*"
                    break
                if nf_std_I < ZERO_EPS:
                    denom = SE_rm
                else:
                    denom = math.sqrt(SE_rm * SE_rm + SE_nf * SE_nf)
                if denom < ZERO_EPS:
                    continue
                tz = dm / denom
                if tz > tz_max:
                    tz_max = tz
                    worst_tz_obs = obs_name

    # -------- Screen and verdict -------------------------------------------
    screen_fail = max_z_mean >= SCREEN_Z_THRESHOLD

    weak = False
    if tz_max is None:
        # No tint data -> fall back to screen-only behavior (historical).
        # This is a weaker metric: screen-only missed example3_fit until
        # 2026-04-10. Surface the downgrade so a reader can tell the
        # verdict came from the fallback rather than the adaptive tz
        # threshold.
        verdict = "FAIL" if screen_fail else "PASS"
        weak = True
        print(
            f"  [weak-verdict] {model}: no tint overlap ({tz_reason}); "
            f"falling back to screen-only (screen={max_z_mean:.2f})",
            file=sys.stderr,
        )
    elif math.isinf(tz_max):
        verdict = "FAIL"
    else:
        verdict = "FAIL" if tz_max >= T else "PASS"

    return {
        "screen_z": max_z_mean,
        "worst_screen": worst_obs_mean,
        "tz_max": tz_max,
        "worst_tz": worst_tz_obs,
        "T_model": T,
        "std_ratio": 1.0 + max_std_ratio,
        "worst_obs_var": worst_obs_var,
        "verdict": verdict,
        "weak": weak,
        "n_compared": n_compared,
    }


# ---------------------------------------------------------------------------
# Markdown output
# ---------------------------------------------------------------------------

MD_HEADER = """\
# RuleMonkey Benchmark Report

**Date:** {date}
**Commit:** `{commit}`
**Reps per model:** {n_reps}
**NFsim reference:** 100-rep ensemble

## Correctness

- **screen**: max |RM_mean - NFsim_mean| / (NFsim_std / sqrt(n_reps)) over all (time, obs) pairs. Fast early-warning metric. Flags values ≥ {screen_thresh} as "suspicious" but does not determine verdict. Historical benchmarks used this as the pass/fail criterion; for models with many rare-event Size_N observables it is statistically unreliable.
- **tz_max**: max over observables of the z-score of the per-rep trapezoidal time integral, computed against precomputed NFsim stats in `ensemble/{{model}}.tint.tsv`. Collapses each 1001-point trajectory to one number per observable per rep, eliminating single-time-point coincidences.
- **T**: per-model verdict threshold = max({tz_floor}, {tz_margin} × tz_p99), where tz_p99 is the 99th percentile of `tz_max` from self-splitting the NFsim replicates at n=10 (see `tests/reference/nfsim/noise_floor.tsv`; provenance and regen recipe in the same directory's `PROVENANCE.md`). Adaptive to each model's intrinsic rare-event noise floor.
- **std_ratio**: max(RM_std / NFsim_std) across observables with nontrivial variance. Diagnostic for variance consistency; not part of verdict.
- **verdict**: PASS if `tz_max < T`, FAIL otherwise. Degenerate-observable mismatches (both stds zero, values differ) fail unconditionally.

## Efficiency

- **nfsim_s**: NFsim mean wall time (100-rep, from reference campaign).
- **rm_s**: RM mean wall time ({n_reps}-rep average).
- **ev/s**: SSA events per wall-second (throughput).
- **sel/fire/obs/upd**: Phase breakdown as % of RM engine time.

## Results

| model | nfsim_s | rm_s | events | ev/s | sel% | fire% | obs% | upd% | screen | tz_max | T | std_ratio | worst_obs | verdict |
|-------|--------:|-----:|-------:|-----:|-----:|------:|-----:|-----:|-------:|-------:|----:|----------:|-----------|---------|
"""

MD_ROW = "| {model} | {nfsim_s} | {rm_s} | {events} | {evs} | {sel} | {fire} | {obs} | {upd} | {screen} | {tz_max} | {T} | {std_ratio} | {worst_obs} | {verdict} |\n"

MD_TIMEOUT_ROW = (
    "| {model} | {nfsim_s} | TIMEOUT | — | — | — | — | — | — | — | — | — | — | — | TIMEOUT |\n"
)
MD_ERROR_ROW = (
    "| {model} | {nfsim_s} | ERROR | — | — | — | — | — | — | — | — | — | — | — | ERROR |\n"
)


def fmt(val, fmt_str=".1f"):
    if val is None:
        return "—"
    return f"{val:{fmt_str}}"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    # Parse args
    args = sys.argv[1:]
    output_path = os.path.join(REPO_ROOT, "build", "benchmark_report.md")

    # --ref-dir DIR re-points all reference paths at a different corpus root.
    # Honored before any other arg parsing so downstream logic sees the new
    # constants. Used by basicmodels.py (and would be used by any future
    # per-corpus driver) to share verdict logic.
    if "--ref-dir" in args:
        idx = args.index("--ref-dir")
        set_ref_dir(args[idx + 1])
        args = args[:idx] + args[idx + 2 :]

    tier = None
    if "--tier" in args:
        idx = args.index("--tier")
        tier = args[idx + 1]
        if tier not in TIER_MODELS:
            print(f"ERROR: --tier must be one of {list(TIER_MODELS.keys())}", file=sys.stderr)
            sys.exit(2)
        args = args[:idx] + args[idx + 2 :]

    # Reps: explicit --reps overrides tier default; tier default overrides
    # global default.
    n_reps = TIER_DEFAULT_REPS.get(tier, DEFAULT_REPS) if tier else DEFAULT_REPS
    if "--reps" in args:
        idx = args.index("--reps")
        n_reps = int(args[idx + 1])
        args = args[:idx] + args[idx + 2 :]
    use_adaptive = "--adaptive" in args
    if use_adaptive:
        args.remove("--adaptive")
    if "--output" in args:
        idx = args.index("--output")
        output_path = args[idx + 1]
        args = args[:idx] + args[idx + 2 :]

    sim_params = load_sim_params()
    nfsim_times = load_nfsim_wall_times()
    model_timeouts = load_timeouts()
    noise_floor = load_noise_floor()
    if not noise_floor:
        print(
            "WARN: no noise_floor.tsv found at "
            f"{NOISE_FLOOR_FILE} — verdict column will fall back to a "
            "flat T=5 threshold (the tsv is a frozen calibration "
            "artifact; see tests/reference/nfsim/PROVENANCE.md).",
            file=sys.stderr,
        )

    # Model list precedence: positional args > tier > all.
    if args:
        models = args
    elif tier and TIER_MODELS[tier] is not None:
        models = [m for m in TIER_MODELS[tier] if m in sim_params]
        missing = [m for m in TIER_MODELS[tier] if m not in sim_params]
        if missing:
            print(
                f"WARN: tier '{tier}' lists models not in sim_params.tsv: {missing}",
                file=sys.stderr,
            )
    else:
        models = sorted(sim_params.keys())

    # Sort fastest-first using NFsim time as proxy, unknown models last
    # (RM-specific slow models will still be slow, but this gets fast ones done first)
    def sort_key(m):
        # Use model_timeouts as best estimate of RM time
        if m in model_timeouts:
            return model_timeouts[m]
        # Otherwise use NFsim time (most models at parity)
        return nfsim_times.get(m, 999999)

    models = sorted(models, key=sort_key)

    # Get git commit
    try:
        commit = subprocess.check_output(
            ["git", "log", "--oneline", "-1"], cwd=REPO_ROOT, text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        commit = "unknown"

    # Check binary exists
    if not os.path.exists(RM_DRIVER):
        print(f"ERROR: {RM_DRIVER} not found. Run: cmake --build build", file=sys.stderr)
        sys.exit(1)

    # Write markdown header
    from datetime import date

    header = MD_HEADER.format(
        date=date.today().isoformat(),
        commit=commit,
        n_reps=n_reps,
        screen_thresh=SCREEN_Z_THRESHOLD,
        tz_floor=VERDICT_TZ_FLOOR,
        tz_margin=VERDICT_TZ_MARGIN,
    )
    with open(output_path, "w") as f:
        f.write(header)

    # Console header
    print(
        f"{'model':<45s} {'nfsim':>7s} {'rm':>7s} {'ev/s':>8s} "
        f"{'screen':>7s} {'tz':>6s} {'T':>5s} {'verdict':<8s}"
    )
    print("-" * 105)

    n_pass = n_fail = n_skip = n_timeout = 0

    for model in models:
        if model not in sim_params:
            n_skip += 1
            continue

        t_end, n_steps, rm_flags = sim_params[model]
        xml = os.path.join(XML_DIR, f"{model}.xml")
        if not os.path.exists(xml):
            n_skip += 1
            continue

        nfsim_s = nfsim_times.get(model)

        # Per-model reps: adaptive mode inflates reps for noisy models
        model_reps = n_reps
        if use_adaptive:
            model_reps = adaptive_reps(noise_floor.get(model), n_reps)

        # Timeout: use model-specific if available, else 10x NFsim time, min 60s
        # Scale timeout by model_reps/10 when adaptive gives more reps
        rep_scale = max(1, model_reps / 10)
        if model in model_timeouts:
            timeout = int(model_timeouts[model] * 2 * rep_scale)
        elif nfsim_s:
            timeout = max(60, int(nfsim_s * 50 * rep_scale))
        else:
            timeout = 3600  # 1 hour default

        # Run reps in parallel (ThreadPoolExecutor; per-rep wall-times are
        # modestly inflated under contention but total model wall scales
        # ~linearly with workers).  Ordered by seed in the result list so
        # downstream stats are deterministic across runs.
        workers = max(1, min(model_reps, MAX_PARALLEL_REPS))
        rep_results = [None] * model_reps
        timed_out = False
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(
                    run_one_rep, xml, t_end, n_steps, 42 + i, timeout=timeout, rm_flags=rm_flags
                ): i
                for i in range(model_reps)
            }
            for fut in as_completed(futures):
                rep_idx = futures[fut]
                result = fut.result()
                if result is None:
                    timed_out = True
                    continue
                rep_results[rep_idx] = result

        all_reps = [r for r in rep_results if r is not None]
        wall_times = [r[3] for r in all_reps]
        timing_info = parse_timing(all_reps[0][2]) if all_reps else None

        if timed_out or len(all_reps) == 0:
            n_timeout += 1
            nf_str = fmt(nfsim_s)
            with open(output_path, "a") as f:
                f.write(MD_TIMEOUT_ROW.format(model=model, nfsim_s=nf_str))
            print(
                f"{model:<45s} {fmt(nfsim_s):>7s} {'T/O':>7s} {'':>8s} {'':>7s} {'':>6s} {'TIMEOUT':<8s}"
            )
            continue

        # Compute results
        rm_s = statistics.mean(wall_times)
        events = timing_info["events"]
        evs = int(events / rm_s) if rm_s > 0 else 0

        # Correctness
        corr = analyze_correctness(model, all_reps, model_reps, nf_row=noise_floor.get(model))

        verdict = corr.get("verdict", "ERROR")
        weak = corr.get("weak", False)
        screen_z = corr.get("screen_z", None)
        tz_val = corr.get("tz_max", None)
        T_val = corr.get("T_model", None)
        std_ratio = corr.get("std_ratio", None)
        worst_tz = corr.get("worst_tz") or corr.get("worst_screen") or ""

        if verdict == "PASS":
            n_pass += 1
        elif verdict == "FAIL":
            n_fail += 1
        else:
            n_skip += 1

        # Annotate verdict with WEAK marker for screen-only fallback.
        # Counters above use the raw PASS/FAIL so tallies are unchanged.
        verdict_display = f"{verdict} (WEAK)" if weak else verdict

        def fmt_tz(v):
            if v is None:
                return "—"
            if isinstance(v, float) and math.isinf(v):
                return "∞"
            return f"{v:.2f}"

        # Write markdown row
        row = MD_ROW.format(
            model=model,
            nfsim_s=fmt(nfsim_s),
            rm_s=fmt(rm_s),
            events=f"{events:,}" if events else "—",
            evs=f"{evs:,}" if evs else "—",
            sel=fmt(timing_info["select_pct"]),
            fire=fmt(timing_info["fire_pct"]),
            obs=fmt(timing_info["obs_pct"]),
            upd=fmt(timing_info["update_pct"]),
            screen=fmt(screen_z, ".2f") if screen_z is not None else "—",
            tz_max=fmt_tz(tz_val),
            T=fmt(T_val, ".2f") if T_val is not None else "—",
            std_ratio=fmt(std_ratio, ".2f") if std_ratio is not None else "—",
            worst_obs=worst_tz[:20] if worst_tz else "—",
            verdict=verdict_display,
        )
        with open(output_path, "a") as f:
            f.write(row)

        # Console output
        reps_str = f"{model_reps}r" if model_reps != n_reps else ""
        print(
            f"{model:<45s} {fmt(nfsim_s):>7s} {fmt(rm_s):>7s} "
            f"{evs:>8,d} {fmt(screen_z, '.2f'):>7s} "
            f"{fmt_tz(tz_val):>6s} {fmt(T_val, '.2f'):>5s} {verdict_display:<8s} {reps_str}"
        )
        sys.stdout.flush()

    # Write summary footer
    footer = f"""
## Summary

| Metric | Count |
|--------|------:|
| PASS | {n_pass} |
| FAIL | {n_fail} |
| TIMEOUT | {n_timeout} |
| SKIP | {n_skip} |
| **Total** | **{n_pass + n_fail + n_timeout + n_skip}** |
"""
    with open(output_path, "a") as f:
        f.write(footer)

    print("-" * 105)
    print(f"PASS: {n_pass}  FAIL: {n_fail}  TIMEOUT: {n_timeout}  SKIP: {n_skip}")
    print(f"\nReport written to: {output_path}")


if __name__ == "__main__":
    main()
