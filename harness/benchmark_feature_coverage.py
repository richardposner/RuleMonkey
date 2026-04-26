#!/usr/bin/env python3
"""Feature-coverage benchmark suite for RuleMonkey.

Tests BNGL features via small, fast models (≤600 molecules, <5s each)
plus a stress tier at larger scale to catch size-dependent bugs.
Designed for fast iteration during development.

Four tiers:
  - base:         One small model per BNGL feature / simulation pattern
  - combinations: Feature interactions that have caused past bugs
  - network-free: Models requiring network-free simulation (NFsim ref only)
  - stress:       Larger-scale, longer-sim-time models exercising incremental
                  bookkeeping (rings, long chains, symmetric complexes,
                  branched aggregates). All network-free.

Reference data:
  - Default (fast): NFsim refs from tests/models/feature_coverage/reference/nfsim/
  - --full mode:    Also reads ODE/SSA refs for cross-checking
  - References are vendored in-tree at tests/models/feature_coverage/reference/.
    If a ref is missing, the script tries to regenerate it via NFSIM_BIN and
    BNG2 (env-overridable). The cleanroom does not vendor those binaries —
    set NFSIM_BIN and BNG2 explicitly when adding a new model. Otherwise
    the missing-ref model is skipped with a warning.

Usage:
  python3 harness/benchmark_feature_coverage.py                  # fast: ~90s
  python3 harness/benchmark_feature_coverage.py --reps 5         # more RM reps
  python3 harness/benchmark_feature_coverage.py --full           # + ODE/SSA refs
  python3 harness/benchmark_feature_coverage.py --tier base      # base only
  python3 harness/benchmark_feature_coverage.py ft_bond_wildcards combo_strict_product_plus

  # Reference regeneration (requires NFSIM_BIN and BNG2 binaries):
  python3 harness/benchmark_feature_coverage.py --force-refs <model>     # regen + run
  python3 harness/benchmark_feature_coverage.py --generate-refs <model>  # regen only, skip RM
"""

import argparse
import glob
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))

# rm_driver: built from this repo. Override via RM_DRIVER env var.
RM_DRIVER = os.environ.get(
    "RM_DRIVER",
    os.path.join(REPO_ROOT, "build", "release", "rm_driver"),
)

# NFsim and BNG2 are external. The cleanroom does not vendor them. They
# are only needed when regenerating references for a model whose ref
# files don't already exist in tests/models/feature_coverage/reference/.
# Defaults point at the sibling nfsim-rm checkout for convenience.
NFSIM = os.environ.get(
    "NFSIM_BIN",
    os.path.expanduser("~/Code/nfsim-rm/build/NFsim"),
)
BNG2 = os.environ.get(
    "BNG2",
    os.path.expanduser("~/Simulations/BioNetGen-2.9.3/BNG2.pl"),
)

SUITE_DIR = os.path.join(REPO_ROOT, "tests", "models", "feature_coverage")
XML_DIR = os.path.join(SUITE_DIR, "xml")
REF_DIR = os.path.join(SUITE_DIR, "reference")
ODE_DIR = os.path.join(REF_DIR, "ode")
SSA_DIR = os.path.join(REF_DIR, "ssa")
NFSIM_DIR = os.path.join(REF_DIR, "nfsim")

REPORT_PATH = os.path.join(REPO_ROOT, "build", "feature_coverage_report.md")

DEFAULT_REPS = 5
DEFAULT_NFSIM_REPS = 20  # enough for stable z-scores with small models
DEFAULT_SSA_REPS = 10
MAX_PARALLEL = 4  # parallel workers for NFsim ref and RM reps

Z_THRESHOLD = 5.0
ZERO_EPS = 1e-9


# ---------------------------------------------------------------------------
# Model discovery and metadata
# ---------------------------------------------------------------------------


def discover_models():
    """Find all .bngl files in the feature_coverage directory."""
    models = []
    for f in sorted(glob.glob(os.path.join(SUITE_DIR, "*.bngl"))):
        models.append(os.path.basename(f)[:-5])
    return models


def model_tier(name):
    if name.startswith("ft_"):
        return "base"
    elif name.startswith("combo_"):
        return "combinations"
    elif name.startswith("nf_"):
        return "network-free"
    elif name.startswith("ss_"):
        return "stress"
    return "unknown"


def is_network_free_only(name):
    return name.startswith("nf_") or name.startswith("ss_")


# Per-tier wall-clock timeouts (seconds).  Stress models run longer
# sims at larger scale, so they need a larger budget than the fast tiers.
def tier_timeout(name):
    if name.startswith("ss_"):
        return 180
    return 30


# Models where NFsim doesn't implement the tested feature correctly.
# For these, we use ODE reference (via BNG2 generate_network) as the
# verdict instead of NFsim.  BNG2's network generation handles
# exclude/include_reactants/products via rule expansion.
NFSIM_UNRELIABLE = {
    "ft_include_reactants",
    "ft_exclude_reactants",
    "ft_exclude_products",
    "combo_exclude_with_complex",
    # NFsim can't resolve function-calling-function; use ODE verdict instead.
    "ft_nested_functions",
    # NFsim (pinned release) silently ignores Fixed species; fix landed
    # upstream in NFsim PR #60 (not in our pinned release).  BNG2 ODE
    # correctly handles Fixed via species_deriv = 0.
    "ft_clamped_species_strict",
}

# Models that declare BNGL features RM refuses by default (Tier-0
# refusal) but whose .bngl additionally provides a workaround path
# that produces a correct trajectory.  For these, the benchmark
# invokes RM with --ignore-unsupported so the comparison against
# NFsim still exercises the non-refused code path.  Verifying that RM
# actually refuses by default is covered by the per-feature test
# fixtures under dev/test_*.xml, not this script.
#
# Currently empty: no corpus or feature-coverage model trips a Tier-0
# trigger (cBNGL: no compartments declared anywhere; eBNGL: only
# Arrhenius rate laws trigger, and nothing in either suite uses them).
TIER0_IGNORE_UNSUPPORTED = set()


def parse_invariants(bngl_path):
    """Extract invariant declarations from BNGL comment block.

    Recognized forms:
      # invariant: nonneg                          (default: always on)
      # invariant: conserved X_total               (matches initial RM value)
      # invariant: conserved X_total = 100         (matches explicit value)
      # invariant: balance R_bound = L_bound       (two observables stay equal)
    """
    inv = {"nonneg": True, "conserved": [], "balance": []}
    with open(bngl_path) as f:
        for line in f:
            m = re.match(r"#\s*invariant\s*:\s*(.+?)\s*$", line.strip())
            if not m:
                continue
            spec = m.group(1).strip()
            if spec == "nonneg":
                inv["nonneg"] = True
            elif spec.startswith("conserved "):
                rest = spec[len("conserved ") :].strip()
                if "=" in rest:
                    name, val = rest.split("=", 1)
                    inv["conserved"].append((name.strip(), float(val.strip())))
                else:
                    inv["conserved"].append((rest, None))
            elif spec.startswith("balance "):
                rest = spec[len("balance ") :].strip()
                if "=" in rest:
                    a, b = rest.split("=", 1)
                    inv["balance"].append((a.strip(), b.strip()))
    return inv


def check_invariants(rm_reps, inv):
    """Validate invariants against RM output. Returns list of violation strings."""
    if not rm_reps:
        return []
    violations = []
    headers = rm_reps[0][0]
    col_idx = {h.strip(): i for i, h in enumerate(headers)}

    for rep_idx, (_h, rows, _e, _w) in enumerate(rm_reps):
        if not rows:
            continue

        if inv.get("nonneg"):
            done = False
            for row in rows:
                for ci in range(1, len(row)):
                    if row[ci] < -ZERO_EPS:
                        violations.append(
                            f"rep{rep_idx}: {headers[ci]}={row[ci]:.3g} < 0 at t={row[0]:.3g}"
                        )
                        done = True
                        break
                if done:
                    break

        for name, expected in inv.get("conserved", []):
            if name not in col_idx:
                violations.append(f"conserved obs '{name}' not found in RM output")
                continue
            ci = col_idx[name]
            init = expected if expected is not None else rows[0][ci]
            tol = max(0.5, 0.01 * abs(init))
            for row in rows:
                if abs(row[ci] - init) > tol:
                    violations.append(
                        f"rep{rep_idx}: {name}={row[ci]:.3g} != expected {init:.3g} "
                        f"(tol {tol:.3g}) at t={row[0]:.3g}"
                    )
                    break

        for a, b in inv.get("balance", []):
            if a not in col_idx or b not in col_idx:
                violations.append(f"balance obs '{a}' or '{b}' not found in RM output")
                continue
            ca, cb = col_idx[a], col_idx[b]
            for row in rows:
                va, vb = row[ca], row[cb]
                tol = max(0.5, 0.01 * max(abs(va), abs(vb)))
                if abs(va - vb) > tol:
                    violations.append(
                        f"rep{rep_idx}: {a}={va:.3g} != {b}={vb:.3g} "
                        f"(tol {tol:.3g}) at t={row[0]:.3g}"
                    )
                    break

    return violations


def parse_model_features(bngl_path):
    """Extract feature declarations from BNGL comment block."""
    features = []
    with open(bngl_path) as f:
        for line in f:
            line = line.strip()
            if not line.startswith("#"):
                break
            if line.startswith("# Feature"):
                features.append(line[2:].strip())
            elif line.startswith("# Tests:"):
                features.append(line[2:].strip())
            elif line.startswith("# Bug pattern:"):
                features.append(line[2:].strip())
    return features


def extract_sim_params(bngl_path):
    """Extract t_end and n_steps from BNGL simulate action."""
    with open(bngl_path) as f:
        text = f.read()
    t_end = None
    n_steps = None
    nfsim_flags = []

    # Find simulate_nf or simulate with method=>"nf"
    for m in re.finditer(r"simulate\w*\s*\(\s*\{([^}]+)\}", text, re.DOTALL):
        args = m.group(1)
        t_match = re.search(r"t_end\s*=>\s*([0-9.eE+\-]+)", args)
        n_match = re.search(r"n_steps\s*=>\s*([0-9]+)", args)
        if t_match:
            t_end = t_match.group(1)
        if n_match:
            n_steps = n_match.group(1)
        # Check for NFsim flags in param=>
        p_match = re.search(r'param\s*=>\s*"([^"]*)"', args)
        if p_match:
            flags_str = p_match.group(1)
            if "-bscb" in flags_str:
                nfsim_flags.append("-bscb")
            utl_match = re.search(r"-utl\s+(\d+)", flags_str)
            if utl_match:
                nfsim_flags.extend(["-utl", utl_match.group(1)])

    return t_end or "100", n_steps or "100", nfsim_flags


# ---------------------------------------------------------------------------
# XML generation
# ---------------------------------------------------------------------------


def generate_xml(model):
    """Generate XML from BNGL via BNG2.pl."""
    os.makedirs(XML_DIR, exist_ok=True)
    xml_path = os.path.join(XML_DIR, f"{model}.xml")
    if os.path.exists(xml_path):
        return xml_path

    bngl_path = os.path.join(SUITE_DIR, f"{model}.bngl")
    with tempfile.TemporaryDirectory() as tmpdir:
        shutil.copy(bngl_path, tmpdir)
        result = subprocess.run(
            ["perl", BNG2, f"{model}.bngl"], cwd=tmpdir, capture_output=True, text=True, timeout=30
        )
        # Find generated XML
        xml_files = glob.glob(os.path.join(tmpdir, "*.xml"))
        if not xml_files:
            print(f"  FAIL: no XML for {model}")
            if result.stderr:
                print(f"  stderr: {result.stderr[:200]}")
            return None
        shutil.copy(xml_files[0], xml_path)
    return xml_path


# ---------------------------------------------------------------------------
# Reference generation
# ---------------------------------------------------------------------------


def generate_ode_reference(model, t_end, n_steps):
    """Generate ODE reference via BNG2.pl generate_network + simulate_ode."""
    os.makedirs(ODE_DIR, exist_ok=True)
    out_path = os.path.join(ODE_DIR, f"{model}.gdat")
    if os.path.exists(out_path):
        return out_path

    bngl_path = os.path.join(SUITE_DIR, f"{model}.bngl")
    with open(bngl_path) as f:
        text = f.read()

    # Replace actions block
    text = re.sub(
        r"begin\s+actions.*?end\s+actions",
        f"begin actions\n"
        f"generate_network({{overwrite=>1}})\n"
        f'simulate({{method=>"ode",t_end=>{t_end},n_steps=>{n_steps}}})\n'
        f"end actions",
        text,
        flags=re.DOTALL,
    )

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_bngl = os.path.join(tmpdir, f"{model}.bngl")
            with open(tmp_bngl, "w") as f:
                f.write(text)
            result = subprocess.run(
                ["perl", BNG2, f"{model}.bngl"],
                cwd=tmpdir,
                capture_output=True,
                text=True,
                timeout=60,
            )
            # Find .gdat file
            gdat_files = glob.glob(os.path.join(tmpdir, "*.gdat"))
            if not gdat_files:
                return None
            # Convert to our TSV format (BNG2 gdat uses spaces)
            h, rows = parse_gdat_file(gdat_files[0])
            write_tsv(out_path, h, rows)
        return out_path
    except subprocess.TimeoutExpired:
        return None


def generate_ssa_reference(model, t_end, n_steps, n_reps=DEFAULT_SSA_REPS):
    """Generate SSA ensemble reference via BNG2.pl."""
    os.makedirs(SSA_DIR, exist_ok=True)
    mean_path = os.path.join(SSA_DIR, f"{model}.mean.tsv")
    std_path = os.path.join(SSA_DIR, f"{model}.std.tsv")
    if os.path.exists(mean_path) and os.path.exists(std_path):
        return mean_path, std_path

    bngl_path = os.path.join(SUITE_DIR, f"{model}.bngl")
    with open(bngl_path) as f:
        text = f.read()

    all_data = []
    headers = None

    for rep in range(n_reps):
        seed = rep + 1
        # Replace actions block with SSA
        rep_text = re.sub(
            r"begin\s+actions.*?end\s+actions",
            f"begin actions\n"
            f"generate_network({{overwrite=>1}})\n"
            f'simulate({{method=>"ssa",t_end=>{t_end},n_steps=>{n_steps},seed=>{seed}}})\n'
            f"end actions",
            text,
            flags=re.DOTALL,
        )

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmp_bngl = os.path.join(tmpdir, f"{model}.bngl")
                with open(tmp_bngl, "w") as f:
                    f.write(rep_text)
                result = subprocess.run(
                    ["perl", BNG2, f"{model}.bngl"],
                    cwd=tmpdir,
                    capture_output=True,
                    text=True,
                    timeout=60,
                )
                gdat_files = glob.glob(os.path.join(tmpdir, "*.gdat"))
                if not gdat_files:
                    continue
                h, rows = parse_gdat_file(gdat_files[0])
                if headers is None:
                    headers = h
                all_data.append(rows)
        except subprocess.TimeoutExpired:
            if rep == 0:
                break  # don't waste time on remaining reps
            continue

    if not all_data or headers is None:
        return None, None

    # Compute mean and std
    n_times = len(all_data[0])
    n_cols = len(all_data[0][0])
    mean_rows = []
    std_rows = []
    for ti in range(n_times):
        mean_row = []
        std_row = []
        for ci in range(n_cols):
            vals = [all_data[ri][ti][ci] for ri in range(len(all_data)) if ti < len(all_data[ri])]
            m = statistics.mean(vals) if vals else 0
            s = statistics.stdev(vals) if len(vals) > 1 else 0
            mean_row.append(m)
            std_row.append(s)
        mean_rows.append(mean_row)
        std_rows.append(std_row)

    write_tsv(mean_path, headers, mean_rows)
    write_tsv(std_path, headers, std_rows)
    return mean_path, std_path


def _run_one_nfsim_rep(model, xml_path, t_end, n_steps, nfsim_flags, seed, timeout=30):
    """Run a single NFsim rep. Returns (headers, rows) or None."""
    with tempfile.TemporaryDirectory() as tmpdir:
        out_gdat = os.path.join(tmpdir, f"{model}_rep{seed}.gdat")
        cmd = [
            NFSIM,
            "-xml",
            xml_path,
            "-sim",
            str(t_end),
            "-oSteps",
            str(n_steps),
            "-seed",
            str(seed),
            "-bscb",
        ]
        for flag in nfsim_flags:
            if flag == "-bscb":
                continue
            cmd.append(flag)
        cmd.extend(["-o", out_gdat])
        try:
            subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
            if not os.path.exists(out_gdat):
                return None
            return parse_gdat_file(out_gdat)
        except subprocess.TimeoutExpired:
            return None


def generate_nfsim_reference(model, t_end, n_steps, nfsim_flags, n_reps=DEFAULT_NFSIM_REPS):
    """Generate NFsim ensemble reference (parallel)."""
    os.makedirs(NFSIM_DIR, exist_ok=True)
    mean_path = os.path.join(NFSIM_DIR, f"{model}.mean.tsv")
    std_path = os.path.join(NFSIM_DIR, f"{model}.std.tsv")
    if os.path.exists(mean_path) and os.path.exists(std_path):
        return mean_path, std_path

    xml_path = os.path.join(XML_DIR, f"{model}.xml")
    if not os.path.exists(xml_path):
        return None, None

    all_data = []
    headers = None

    nf_timeout = tier_timeout(model)
    with ThreadPoolExecutor(max_workers=MAX_PARALLEL) as pool:
        futures = {
            pool.submit(
                _run_one_nfsim_rep, model, xml_path, t_end, n_steps, nfsim_flags, seed, nf_timeout
            ): seed
            for seed in range(1, n_reps + 1)
        }
        for future in as_completed(futures):
            result = future.result()
            if result is not None:
                h, rows = result
                if headers is None:
                    headers = h
                all_data.append(rows)

    if not all_data or headers is None:
        return None, None

    # Compute mean and std (guard against empty rep data)
    if not all_data[0] or not all_data[0][0]:
        return None, None
    n_times = len(all_data[0])
    n_cols = len(all_data[0][0])
    mean_rows = []
    std_rows = []
    for ti in range(n_times):
        mean_row = []
        std_row = []
        for ci in range(n_cols):
            vals = [all_data[ri][ti][ci] for ri in range(len(all_data)) if ti < len(all_data[ri])]
            m = statistics.mean(vals) if vals else 0
            s = statistics.stdev(vals) if len(vals) > 1 else 0
            mean_row.append(m)
            std_row.append(s)
        mean_rows.append(mean_row)
        std_rows.append(std_row)

    write_tsv(mean_path, headers, mean_rows)
    write_tsv(std_path, headers, std_rows)
    return mean_path, std_path


# ---------------------------------------------------------------------------
# Parsing & I/O
# ---------------------------------------------------------------------------


def parse_gdat(text):
    """Parse gdat text output (from stdout)."""
    lines = [l.strip() for l in text.strip().split("\n") if l.strip()]
    # Find header line (first line starting with # or containing 'time')
    header_line = lines[0]
    headers = header_line.lstrip("#").strip().split()
    # Clean up headers: some NFsim outputs use tabs, some spaces
    if len(headers) <= 1:
        headers = header_line.lstrip("#").strip().split("\t")
    rows = []
    for line in lines[1:]:
        if line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) <= 1:
            parts = line.split("\t")
        try:
            rows.append([float(v) for v in parts])
        except ValueError:
            continue
    return headers, rows


def parse_gdat_file(path):
    """Parse gdat file."""
    with open(path) as f:
        # BNG2.pl gdat uses whitespace-separated columns with # header
        lines = [l.strip() for l in f if l.strip()]
    header_line = lines[0]
    headers = header_line.lstrip("#").strip().split()
    if len(headers) <= 1:
        headers = header_line.lstrip("#").strip().split("\t")
    rows = []
    for line in lines[1:]:
        if line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) <= 1:
            parts = line.split("\t")
        try:
            rows.append([float(v) for v in parts])
        except ValueError:
            continue
    return headers, rows


def read_tsv(path):
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip()]
    headers = lines[0].split("\t")
    rows = [[float(v) for v in line.split("\t")] for line in lines[1:]]
    return headers, rows


def write_tsv(path, headers, rows):
    with open(path, "w") as f:
        f.write("\t".join(headers) + "\n")
        for row in rows:
            f.write("\t".join(f"{v:.6g}" for v in row) + "\n")


def parse_timing(stderr_text):
    """Extract timing info from RM stderr output."""
    info = {"events": 0, "total_s": 0}
    m = re.search(r"events=(\d+)", stderr_text)
    if m:
        info["events"] = int(m.group(1))
    m = re.search(r"total=([\d.]+)s", stderr_text)
    if m:
        info["total_s"] = float(m.group(1))
    return info


# ---------------------------------------------------------------------------
# RM execution
# ---------------------------------------------------------------------------


def run_rm_rep(xml_path, t_end, n_steps, seed, rm_flags=None, timeout=30):
    """Run a single RM rep."""
    cmd = [RM_DRIVER, xml_path, str(t_end), str(n_steps), str(seed)]
    if rm_flags:
        cmd.extend(rm_flags)
    t0 = time.monotonic()
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        wall_s = time.monotonic() - t0
        if result.returncode != 0 or not result.stdout.strip():
            return None
        headers, rows = parse_gdat(result.stdout)
        return (headers, rows, result.stderr, wall_s)
    except (subprocess.TimeoutExpired, FileNotFoundError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Correctness comparison
# ---------------------------------------------------------------------------


def compare_rm_vs_reference(model, rm_reps, ref_mean_path, ref_std_path, ref_label="ref"):
    """Compare RM ensemble against a reference (ODE, SSA, or NFsim).

    Returns dict with max_z, worst_obs, tz_max, verdict.
    """
    if not os.path.exists(ref_mean_path):
        return {"verdict": "NO_REF", "label": ref_label}
    if ref_std_path is None or not os.path.exists(ref_std_path):
        # ODE has no std — use mean-only comparison
        ref_has_std = False
    else:
        ref_has_std = True

    ref_headers, ref_mean = read_tsv(ref_mean_path)
    if ref_has_std:
        _, ref_std = read_tsv(ref_std_path)
    else:
        ref_std = None

    rm_headers = rm_reps[0][0]
    n_reps = len(rm_reps)

    # Column mapping: RM header -> ref column index
    col_map = {}
    for ri, rh in enumerate(rm_headers):
        for ni, nh in enumerate(ref_headers):
            if rh.strip() == nh.strip():
                col_map[ri] = ni
                break

    # Ref time -> row index
    ref_time_idx = {round(row[0], 6): i for i, row in enumerate(ref_mean)}

    # Collect per-timepoint values
    max_z = 0
    worst_obs = ""
    n_compared = 0
    max_abs_diff = 0
    worst_abs_obs = ""

    # Get shared time points
    time_points = []
    for rm_row in rm_reps[0][1]:
        rm_t = round(rm_row[0], 6)
        if rm_t in ref_time_idx:
            time_points.append(rm_t)

    # Collect RM values per (time, obs)
    obs_data = {}
    for rm_ci, ref_ci in col_map.items():
        if rm_ci == 0 or ref_ci == 0:
            continue
        obs_name = rm_headers[rm_ci]
        obs_data[obs_name] = {tp: [] for tp in time_points}

    for rep_h, rep_rows, _, _ in rm_reps:
        for rm_row in rep_rows:
            rm_t = round(rm_row[0], 6)
            if rm_t not in ref_time_idx:
                continue
            for rm_ci, ref_ci in col_map.items():
                if rm_ci == 0 or ref_ci == 0:
                    continue
                obs_name = rm_headers[rm_ci]
                obs_data[obs_name][rm_t].append(rm_row[rm_ci])

    for obs_name, tdata in obs_data.items():
        for tp, rm_vals in tdata.items():
            if not rm_vals:
                continue
            n_compared += 1
            ref_idx = ref_time_idx[tp]

            ref_ci = None
            for rm_ci, rc in col_map.items():
                if rm_headers[rm_ci] == obs_name:
                    ref_ci = rc
                    break
            if ref_ci is None:
                continue

            ref_m = ref_mean[ref_idx][ref_ci]
            rm_m = statistics.mean(rm_vals)

            # Absolute difference (for ODE comparison where std may not exist)
            abs_diff = abs(rm_m - ref_m)
            if abs_diff > max_abs_diff:
                max_abs_diff = abs_diff
                worst_abs_obs = obs_name

            # Z-score (only when we have std)
            if ref_has_std and ref_std is not None:
                ref_s = ref_std[ref_idx][ref_ci]
                if ref_s > ZERO_EPS:
                    z = abs(rm_m - ref_m) / (ref_s / math.sqrt(n_reps))
                elif abs(rm_m - ref_m) < max(1.0, 0.05 * (abs(ref_m) + 1.0)):
                    z = 0.0  # both near-constant and close enough
                else:
                    z = 100.0  # std=0 but values differ significantly
                if z > max_z:
                    max_z = z
                    worst_obs = obs_name

    # Compute time-integrated z-score (tz_max)
    tz_max = 0
    worst_tz_obs = ""
    tp_sorted = sorted(time_points)

    if len(tp_sorted) >= 2 and ref_has_std:
        for rm_ci, ref_ci in col_map.items():
            if rm_ci == 0 or ref_ci == 0:
                continue
            obs_name = rm_headers[rm_ci]

            # Per-rep trapezoidal integrals
            per_rep_I = []
            for rep_h, rep_rows, _, _ in rm_reps:
                t2v = {round(r[0], 6): r for r in rep_rows}
                vals = []
                ok = True
                for tp in tp_sorted:
                    if tp not in t2v:
                        ok = False
                        break
                    vals.append(t2v[tp][rm_ci])
                if not ok:
                    continue
                I = sum(
                    0.5 * (vals[i] + vals[i + 1]) * (tp_sorted[i + 1] - tp_sorted[i])
                    for i in range(len(tp_sorted) - 1)
                )
                per_rep_I.append(I)

            if not per_rep_I:
                continue

            # Ref integral (trapezoidal on ref mean)
            ref_vals = []
            for tp in tp_sorted:
                ref_vals.append(ref_mean[ref_time_idx[tp]][ref_ci])
            I_ref = sum(
                0.5 * (ref_vals[i] + ref_vals[i + 1]) * (tp_sorted[i + 1] - tp_sorted[i])
                for i in range(len(tp_sorted) - 1)
            )

            rm_I_mean = statistics.mean(per_rep_I)
            rm_I_se = (
                statistics.stdev(per_rep_I) / math.sqrt(len(per_rep_I)) if len(per_rep_I) > 1 else 0
            )

            denom = rm_I_se if rm_I_se > ZERO_EPS else 1.0
            tz = abs(rm_I_mean - I_ref) / denom if denom > ZERO_EPS else 0
            if tz > tz_max:
                tz_max = tz
                worst_tz_obs = obs_name

    # Verdict: tz_max is the primary metric (robust); max_z is diagnostic only
    # (max_z is noisy with few reps due to multiple-comparison inflation)
    verdict = "PASS" if tz_max < Z_THRESHOLD else "FAIL"
    # For ODE (no std), use relative error with generous threshold
    # (stochastic RM will naturally differ from deterministic ODE by ~sqrt(N))
    if not ref_has_std:
        max_ref_val = (
            max(
                abs(ref_mean[i][j])
                for i in range(len(ref_mean))
                for j in range(1, len(ref_mean[0]))
            )
            if ref_mean
            else 1.0
        )
        rel_err = max_abs_diff / (max_ref_val + 1.0)
        verdict = "PASS" if rel_err < 0.25 else "FAIL"

    return {
        "label": ref_label,
        "max_z": max_z,
        "worst_obs": worst_obs,
        "tz_max": tz_max,
        "worst_tz_obs": worst_tz_obs,
        "max_abs_diff": max_abs_diff,
        "n_compared": n_compared,
        "verdict": verdict,
    }


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------


def write_report(results, report_path):
    """Write markdown benchmark report."""
    with open(report_path, "w") as f:
        f.write("# Feature Coverage Benchmark Report\n\n")
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        # Summary
        n_pass = sum(1 for r in results if r["verdict"] == "PASS")
        n_fail = sum(1 for r in results if r["verdict"] == "FAIL")
        n_skip = sum(1 for r in results if r["verdict"] in ("SKIP", "ERROR", "NO_REF"))
        f.write(f"**Summary: {n_pass} PASS / {n_fail} FAIL / {n_skip} SKIP**\n\n")

        # Coverage matrix
        f.write("## Feature Coverage\n\n")
        f.write("| Model | Tier | Features | RM vs NFsim | RM vs ODE | Verdict |\n")
        f.write("|-------|------|----------|-------------|-----------|--------|\n")

        for r in results:
            tier = model_tier(r["model"])
            features = r.get("features", [""])
            feat_str = features[0][:50] if features else ""
            nfsim_str = r.get("nfsim_verdict", "-")
            ode_str = r.get("ode_verdict", "-")
            verdict = r["verdict"]
            mark = "PASS" if verdict == "PASS" else f"**{verdict}**"
            f.write(f"| {r['model']} | {tier} | {feat_str} | {nfsim_str} | {ode_str} | {mark} |\n")

        # Detailed results
        f.write("\n## Detailed Results\n\n")
        for r in results:
            if r["verdict"] in ("SKIP", "ERROR"):
                f.write(f"### {r['model']} — {r['verdict']}\n")
                if "error" in r:
                    f.write(f"  Error: {r['error']}\n\n")
                continue

            f.write(f"### {r['model']}\n")
            f.write(f"- Tier: {model_tier(r['model'])}\n")
            f.write(f"- RM reps: {r.get('n_reps', '?')}, wall time: {r.get('rm_wall_s', 0):.3f}s\n")

            for comp in r.get("comparisons", []):
                label = comp["label"]
                if comp["verdict"] == "NO_REF":
                    f.write(f"- vs {label}: no reference data\n")
                else:
                    f.write(f"- vs {label}: max_z={comp['max_z']:.2f}")
                    if comp.get("worst_obs"):
                        f.write(f" ({comp['worst_obs']})")
                    f.write(f", tz_max={comp.get('tz_max', 0):.2f}")
                    f.write(f" — **{comp['verdict']}**\n")

            inv_v = r.get("invariant_violations", [])
            if inv_v:
                f.write(f"- invariant violations: {len(inv_v)}\n")
                for v in inv_v[:5]:
                    f.write(f"  - {v}\n")
                if len(inv_v) > 5:
                    f.write(f"  - ... {len(inv_v) - 5} more\n")

            f.write(f"- **Overall: {r['verdict']}**\n\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def run_benchmark(
    models,
    n_reps,
    generate_refs_only=False,
    nfsim_reps=DEFAULT_NFSIM_REPS,
    ssa_reps=DEFAULT_SSA_REPS,
    full_mode=False,
):
    """Run the feature-coverage benchmark."""
    results = []

    for model in models:
        print(f"\n{'=' * 60}")
        print(f"  {model} ({model_tier(model)})")
        print(f"{'=' * 60}")

        bngl_path = os.path.join(SUITE_DIR, f"{model}.bngl")
        if not os.path.exists(bngl_path):
            print(f"  SKIP: {bngl_path} not found")
            results.append({"model": model, "verdict": "SKIP", "error": "BNGL not found"})
            continue

        features = parse_model_features(bngl_path)
        t_end, n_steps, nfsim_flags = extract_sim_params(bngl_path)

        # 1. Generate XML
        print("  XML...", end=" ", flush=True)
        xml_path = generate_xml(model)
        if not xml_path:
            results.append(
                {
                    "model": model,
                    "verdict": "ERROR",
                    "error": "XML gen failed",
                    "features": features,
                }
            )
            continue
        print("OK", end="", flush=True)

        # 2. Generate references
        nf_only = is_network_free_only(model)
        use_ode_verdict = model in NFSIM_UNRELIABLE
        rm_flags = []
        if "-bscb" in nfsim_flags:
            rm_flags.append("-bscb")
        if model in TIER0_IGNORE_UNSUPPORTED:  # noqa: currently empty
            rm_flags.append("--ignore-unsupported")

        ode_path = None
        ssa_mean = ssa_std = None

        if (full_mode and not nf_only) or use_ode_verdict:
            print(", ODE...", end="", flush=True)
            ode_path = generate_ode_reference(model, t_end, n_steps)
            print("OK" if ode_path else "skip", end="", flush=True)

            if full_mode and not nf_only:
                print(f", SSA({ssa_reps})...", end="", flush=True)
                ssa_mean, ssa_std = generate_ssa_reference(model, t_end, n_steps, ssa_reps)
                print("OK" if ssa_mean else "skip", end="", flush=True)

        print(f", NFsim({nfsim_reps})...", end="", flush=True)
        nf_mean, nf_std = generate_nfsim_reference(model, t_end, n_steps, nfsim_flags, nfsim_reps)
        print("OK" if nf_mean else "skip", flush=True)

        if generate_refs_only:
            results.append({"model": model, "verdict": "REF_ONLY", "features": features})
            continue

        # 3. Run RM (parallel)
        rm_reps = []
        rm_timeout = tier_timeout(model)
        t0 = time.monotonic()
        with ThreadPoolExecutor(max_workers=MAX_PARALLEL) as pool:
            futures = {
                pool.submit(run_rm_rep, xml_path, t_end, n_steps, 42 + i, rm_flags, rm_timeout): i
                for i in range(n_reps)
            }
            for future in as_completed(futures):
                rep = future.result()
                if rep:
                    rm_reps.append(rep)
        rm_wall = time.monotonic() - t0
        print(f"  RM: {len(rm_reps)}/{n_reps} OK, {rm_wall:.2f}s", end="", flush=True)

        if not rm_reps:
            print()
            results.append(
                {
                    "model": model,
                    "verdict": "ERROR",
                    "error": "all RM reps failed",
                    "features": features,
                }
            )
            continue

        # 4. Compare
        comparisons = []

        # vs NFsim (primary verdict)
        if nf_mean and nf_std:
            comp = compare_rm_vs_reference(model, rm_reps, nf_mean, nf_std, "NFsim")
            comparisons.append(comp)
            print(f", tz={comp.get('tz_max', 0):.2f}->{comp['verdict']}", end="", flush=True)

        if (full_mode and not nf_only) or use_ode_verdict:
            # vs ODE (verdict for NFSIM_UNRELIABLE models, informational otherwise)
            if ode_path:
                comp_ode = compare_rm_vs_reference(model, rm_reps, ode_path, None, "ODE")
                comparisons.append(comp_ode)
                if use_ode_verdict:
                    print(
                        f", ODE_tz={comp_ode.get('tz_max', 0):.2f}->{comp_ode['verdict']}",
                        end="",
                        flush=True,
                    )

            # vs SSA (informational)
            if full_mode and not nf_only and ssa_mean and ssa_std:
                comp_ssa = compare_rm_vs_reference(model, rm_reps, ssa_mean, ssa_std, "SSA")
                comparisons.append(comp_ssa)

        # Overall verdict: NFsim is normally the sole verdict reference.
        # For models in NFSIM_UNRELIABLE (where NFsim ignores the tested
        # feature), ODE reference is used as the verdict instead.
        nfsim_comp = next((c for c in comparisons if c["label"] == "NFsim"), None)
        ode_comp = next((c for c in comparisons if c["label"] == "ODE"), None)

        if use_ode_verdict and ode_comp and ode_comp["verdict"] != "NO_REF":
            overall = ode_comp["verdict"]
        elif nfsim_comp and nfsim_comp["verdict"] != "NO_REF":
            overall = nfsim_comp["verdict"]
        else:
            # No NFsim ref: fall back to any available
            verdicts = [c["verdict"] for c in comparisons if c["verdict"] != "NO_REF"]
            overall = (
                "FAIL" if any(v == "FAIL" for v in verdicts) else ("PASS" if verdicts else "NO_REF")
            )

        # 5. Invariant checks (non-negativity always; conservation/balance per model)
        invariants = parse_invariants(bngl_path)
        inv_violations = check_invariants(rm_reps, invariants)
        if inv_violations:
            print(f", INV:{len(inv_violations)}", end="", flush=True)
            overall = "FAIL"

        nfsim_v = next((c["verdict"] for c in comparisons if c["label"] == "NFsim"), "-")
        ode_v = next((c["verdict"] for c in comparisons if c["label"] == "ODE"), "-")

        results.append(
            {
                "model": model,
                "verdict": overall,
                "features": features,
                "comparisons": comparisons,
                "n_reps": len(rm_reps),
                "rm_wall_s": rm_wall,
                "nfsim_verdict": nfsim_v,
                "ode_verdict": ode_v,
                "invariant_violations": inv_violations,
            }
        )

        # Parse RM timing from last rep
        timing = parse_timing(rm_reps[-1][2])
        if timing["events"] > 0:
            print(f", {timing['events']}ev", end="", flush=True)
        print(flush=True)

    return results


def main():
    parser = argparse.ArgumentParser(description="Feature coverage benchmark")
    parser.add_argument("models", nargs="*", help="Specific model names (omit for all)")
    parser.add_argument("--reps", type=int, default=DEFAULT_REPS, help="RM replicate count")
    parser.add_argument(
        "--tier",
        choices=["base", "combinations", "network-free", "stress", "all"],
        default="all",
        help="Run only models from this tier",
    )
    parser.add_argument(
        "--full", action="store_true", help="Full mode: also generate ODE/SSA references (slower)"
    )
    parser.add_argument(
        "--generate-refs",
        action="store_true",
        help="Generate reference data only (no RM comparison)",
    )
    parser.add_argument(
        "--nfsim-reps", type=int, default=DEFAULT_NFSIM_REPS, help="NFsim reference replicate count"
    )
    parser.add_argument(
        "--ssa-reps", type=int, default=DEFAULT_SSA_REPS, help="SSA reference replicate count"
    )
    parser.add_argument("--output", default=REPORT_PATH, help="Output report path")
    parser.add_argument(
        "--force-refs", action="store_true", help="Regenerate reference data even if it exists"
    )
    args = parser.parse_args()

    # Discover models
    all_models = discover_models()
    if args.models:
        models = [m for m in args.models if m in all_models]
        if not models:
            print(f"No matching models found. Available: {', '.join(all_models[:10])}...")
            sys.exit(1)
    elif args.tier != "all":
        tier_prefix = {
            "base": "ft_",
            "combinations": "combo_",
            "network-free": "nf_",
            "stress": "ss_",
        }
        prefix = tier_prefix.get(args.tier, "")
        models = [m for m in all_models if m.startswith(prefix)]
    else:
        models = all_models

    mode = "full (NFsim+ODE+SSA)" if args.full else "fast (NFsim only)"
    print(f"Feature Coverage Benchmark — {mode}")
    print(f"  {len(models)} models, {args.reps} RM reps, {args.nfsim_reps} NFsim ref reps")

    if args.force_refs:
        dirs_to_clear = [NFSIM_DIR]
        if args.full:
            dirs_to_clear.extend([ODE_DIR, SSA_DIR])
        # Only clear refs for the selected models. Full-suite regeneration
        # happens naturally when `models` covers everything.
        for d in dirs_to_clear:
            if not os.path.isdir(d):
                continue
            for m in models:
                for path in glob.glob(os.path.join(d, f"{m}.*")):
                    os.remove(path)

    # Ensure directories
    for d in [XML_DIR, ODE_DIR, SSA_DIR, NFSIM_DIR]:
        os.makedirs(d, exist_ok=True)

    results = run_benchmark(
        models,
        args.reps,
        generate_refs_only=args.generate_refs,
        nfsim_reps=args.nfsim_reps,
        ssa_reps=args.ssa_reps,
        full_mode=args.full,
    )

    # Write report
    write_report(results, args.output)
    print(f"\nReport written to {args.output}")

    # Summary
    n_pass = sum(1 for r in results if r["verdict"] == "PASS")
    n_fail = sum(1 for r in results if r["verdict"] == "FAIL")
    n_other = len(results) - n_pass - n_fail
    print(f"\n{'=' * 60}")
    print(f"  TOTAL: {n_pass} PASS / {n_fail} FAIL / {n_other} other")
    print(f"{'=' * 60}")

    sys.exit(1 if n_fail > 0 else 0)


if __name__ == "__main__":
    main()
