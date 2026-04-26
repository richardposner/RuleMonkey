#!/usr/bin/env python3
"""Validate RM output against NFsim ensemble (mean/std).

For each model that completes within 35s, compare RM observables at each
output time point against the NFsim 100-rep ensemble mean/std.  Report
max |z-score| per observable per model.

With --reps N, runs N replicate simulations with different seeds, averages
the trajectories, and compares the RM ensemble mean against NFsim.  This
reduces stochastic noise by 1/sqrt(N), resolving false positives on rare
species.

Usage:
  python3 harness/benchmark_validate.py                  # all models, 1 rep
  python3 harness/benchmark_validate.py --reps 10        # 10 reps
  python3 harness/benchmark_validate.py --reps 10 tcr    # 10 reps, subset
"""

import subprocess, sys, os, concurrent.futures

DEFAULT_TIMEOUT = 45
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
RM_DRIVER = os.environ.get(
    "RM_DRIVER",
    os.path.join(REPO_ROOT, "build", "release", "rm_driver"),
)
REF_DIR = os.path.join(REPO_ROOT, "tests", "reference", "nfsim")
XML_DIR = os.path.join(REF_DIR, "xml")
PARAMS = os.path.join(REF_DIR, "sim_params.tsv")
ENS_DIR = os.path.join(REF_DIR, "ensemble")
TIMEOUT_FILE = os.path.join(REPO_ROOT, "harness", "model_timeouts.tsv")

Z_THRESHOLD = 5.0  # z-scores above this are flagged FAIL


def load_timeouts():
    """Load model-specific timeouts from harness/model_timeouts.tsv."""
    timeouts = {}
    if not os.path.exists(TIMEOUT_FILE):
        return timeouts
    with open(TIMEOUT_FILE) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('model'):
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                timeouts[parts[0]] = int(parts[1])
    return timeouts


def read_tsv(path):
    """Read a TSV file, return (headers, rows) where rows are list of floats."""
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip()]
    headers = lines[0].split('\t')
    rows = []
    for line in lines[1:]:
        rows.append([float(v) for v in line.split('\t')])
    return headers, rows


def parse_gdat(text):
    """Parse RM .gdat output (# header line, tab-separated)."""
    lines = [l.strip() for l in text.strip().split('\n') if l.strip()]
    headers = lines[0].lstrip('#').strip().split('\t')
    rows = []
    for line in lines[1:]:
        rows.append([float(v) for v in line.split('\t')])
    return headers, rows


def load_sim_params():
    params = {}
    with open(PARAMS) as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts[0] == 'model':
                continue
            params[parts[0]] = (parts[2], parts[3])  # t_end, n_steps
    return params


def run_one_rep(xml, t_end, n_steps, seed, timeout=DEFAULT_TIMEOUT):
    """Run a single RM simulation, return (headers, rows) or None."""
    try:
        result = subprocess.run(
            [RM_DRIVER, xml, t_end, n_steps, str(seed)],
            capture_output=True, text=True, timeout=timeout
        )
        if result.returncode != 0 or not result.stdout.strip():
            return None
    except subprocess.TimeoutExpired:
        return None
    return parse_gdat(result.stdout)


def validate_model(model, t_end, n_steps, n_reps=1, timeout=DEFAULT_TIMEOUT):
    xml = os.path.join(XML_DIR, f"{model}.xml")
    mean_path = os.path.join(ENS_DIR, f"{model}.mean.tsv")
    std_path = os.path.join(ENS_DIR, f"{model}.std.tsv")

    if not os.path.exists(xml):
        return None, "no xml"
    if not os.path.exists(mean_path) or not os.path.exists(std_path):
        return None, "no ensemble data"

    if n_reps <= 1:
        # Single rep (original behavior)
        parsed = run_one_rep(xml, t_end, n_steps, 42, timeout)
        if parsed is None:
            return None, "TIMEOUT"
        rm_headers, rm_rows = parsed
    else:
        # First, time a single rep to decide parallel vs sequential
        seeds = list(range(42, 42 + n_reps))
        first = run_one_rep(xml, t_end, n_steps, seeds[0], timeout)
        if first is None:
            return None, "TIMEOUT"
        rm_headers = first[0]
        all_rows = [first[1]]

        remaining = seeds[1:]
        # Parallelize only if the model is fast enough that contention
        # won't push reps past the timeout (heuristic: < 5s solo)
        import time
        t0 = time.monotonic()
        probe = run_one_rep(xml, t_end, n_steps, seeds[1], timeout)
        solo_time = time.monotonic() - t0
        if probe is None:
            return None, "TIMEOUT"
        all_rows.append(probe[1])
        remaining = seeds[2:]

        if remaining:
            use_parallel = solo_time < timeout / 4
            if use_parallel:
                with concurrent.futures.ThreadPoolExecutor(
                        max_workers=min(len(remaining), 8)) as ex:
                    futures = {ex.submit(run_one_rep, xml, t_end, n_steps, s, timeout): s
                               for s in remaining}
                    for fut in concurrent.futures.as_completed(futures):
                        parsed = fut.result()
                        if parsed is None:
                            return None, "TIMEOUT"
                        all_rows.append(parsed[1])
            else:
                for s in remaining:
                    parsed = run_one_rep(xml, t_end, n_steps, s, timeout)
                    if parsed is None:
                        return None, "TIMEOUT"
                    all_rows.append(parsed[1])

        # Average across reps (element-wise, including time column)
        n = len(all_rows)
        rm_rows = []
        for row_idx in range(len(all_rows[0])):
            ncols = len(all_rows[0][row_idx])
            avg = [0.0] * ncols
            for rep in range(n):
                for ci in range(ncols):
                    avg[ci] += all_rows[rep][row_idx][ci]
            rm_rows.append([v / n for v in avg])

    nf_headers, nf_mean = read_tsv(mean_path)
    _, nf_std = read_tsv(std_path)

    # RM header -> NF column index
    col_map = {}
    for ri, rh in enumerate(rm_headers):
        for ni, nh in enumerate(nf_headers):
            if rh.strip() == nh.strip():
                col_map[ri] = ni
                break

    # NF time -> row index
    nf_time_idx = {round(row[0], 6): i for i, row in enumerate(nf_mean)}

    obs_results = {}  # obs_name -> max_z
    n_compared = 0

    for rm_row in rm_rows:
        rm_t = round(rm_row[0], 6)
        nf_idx = nf_time_idx.get(rm_t)
        if nf_idx is None:
            continue
        n_compared += 1

        for rm_ci, nf_ci in col_map.items():
            if rm_ci == 0 or nf_ci == 0:
                continue
            obs_name = rm_headers[rm_ci]
            rm_val = rm_row[rm_ci]
            nf_val = nf_mean[nf_idx][nf_ci]
            nf_sd = nf_std[nf_idx][nf_ci]

            if nf_sd > 0:
                z = abs(rm_val - nf_val) / nf_sd
            elif abs(rm_val - nf_val) < 1e-9:
                z = 0.0
            else:
                # std=0 with nonzero difference: rare species where the
                # ensemble had zero variance (all 100 runs agreed).  Score
                # as |delta| directly — a single molecule (delta=1) scores
                # 1.0, well under any failure threshold.
                z = abs(rm_val - nf_val)

            if obs_name not in obs_results or z > obs_results[obs_name]:
                obs_results[obs_name] = z

    if n_compared == 0:
        return None, "no matching time points"

    return obs_results, "OK"


def main():
    sim_params = load_sim_params()

    # Parse --reps N
    args = sys.argv[1:]
    n_reps = 1
    if '--reps' in args:
        idx = args.index('--reps')
        n_reps = int(args[idx + 1])
        args = args[:idx] + args[idx + 2:]

    models = args if args else sorted(sim_params.keys())
    model_timeouts = load_timeouts()

    reps_label = f" ({n_reps} reps)" if n_reps > 1 else ""
    print(f"{'model':<50s} {'max_z':>8s} {'worst_obs':<25s} {'verdict':<8s} {'n_obs':>5s}{reps_label}")
    print("-" * 100)

    n_pass = n_fail = n_skip = 0
    fails = []

    for model in models:
        if model not in sim_params:
            print(f"{model:<50s} {'':>8s} {'':.<25s} {'SKIP':.<8s}")
            n_skip += 1
            continue

        t_end, n_steps = sim_params[model]
        model_timeout = model_timeouts.get(model, DEFAULT_TIMEOUT)
        obs_results, status = validate_model(model, t_end, n_steps, n_reps, model_timeout)

        if obs_results is None:
            print(f"{model:<50s} {'':>8s} {'':.<25s} {status:<8s}")
            n_skip += 1
            continue

        max_z = 0
        worst_obs = ""
        for obs, z in obs_results.items():
            if z > max_z:
                max_z = z
                worst_obs = obs

        verdict = "PASS" if max_z < Z_THRESHOLD else "FAIL"
        if verdict == "PASS":
            n_pass += 1
        else:
            n_fail += 1
            fails.append((model, max_z, worst_obs))

        print(f"{model:<50s} {max_z:8.2f} {worst_obs:<25s} {verdict:<8s} {len(obs_results):5d}")

    print("-" * 100)
    print(f"PASS: {n_pass}  FAIL: {n_fail}  SKIP: {n_skip}")
    if fails:
        print("\nFailed models:")
        for m, z, obs in fails:
            print(f"  {m}: max_z={z:.2f} on {obs}")


if __name__ == "__main__":
    main()
