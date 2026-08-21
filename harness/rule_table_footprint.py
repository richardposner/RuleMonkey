#!/usr/bin/env python3
"""Per-rule side-table residency across the three corpora (issue #71).

Every per-molecule table a rule owns is sized by the molecule ARENA and not
by the population the rule can actually see, so what a model holds in them
is a (rule count x arena) product.  #72 narrowed the row -- 80 bytes to 12,
with the wide half allocated only where it is read -- and left the sizing
alone, on the grounds that the case for changing it is residency on real
models rather than on a synthetic scale-up, and that nobody had measured
what real models hold.

This is that measurement.  For each model in `feature_coverage`, `corpus`
and `nfsim_basicmodels` it runs one replicate at the suite's canonical
(t_end, n_steps), and records:

  * `rules` / `tabled`  -- rules in the expanded model, and how many hold
                           a per-molecule table at all;
  * `arena`             -- `pool.molecule_count()` at the end of the run,
                           which is what the tables are sized by;
  * `table_bytes`       -- what those tables hold, split into per-molecule
                           rows and Fenwick samplers;
  * `reach_bytes`       -- what the same tables would hold if each were cut
                           to the highest live molecule id its rule can
                           index.  That is the sizing #71 names as the
                           obvious shortcut and warns does not generalize;
                           the gap against `table_bytes` is what it would
                           buy, per model, measured rather than argued;
  * `peak_rss`          -- peak resident set of the rm_driver process;
  * `share`             -- table_bytes / peak_rss, the number the sizing
                           question turns on.

All but the last two come from the engine's `[RM tables]` stderr block,
which `RM_PRINT_TABLES=1` turns on (see RuleTableBytes in engine.cpp for
what is counted, and why counting rows written rather than rows allocated
makes this a floor on the residency).  The peak RSS is the child process's
own `ru_maxrss`, taken from `os.wait4` so it is that run's number and not a
high-water mark over every child this script has spawned.

Usage:
  python3 harness/rule_table_footprint.py                     # all 189 models
  python3 harness/rule_table_footprint.py --suite corpus      # one suite
  python3 harness/rule_table_footprint.py --jobs 4            # in parallel
  python3 harness/rule_table_footprint.py AN BLBR r07         # named models
  python3 harness/rule_table_footprint.py --output report.md

Serial by default: the RSS of a run is not affected by what else is
running, but its wall time is, and the report carries both.  The full
sweep is roughly seven minutes serial.

POSIX only, unlike the parity harnesses: `os.wait4` is where the child's
own peak RSS comes from, and there is no Windows equivalent that does not
mean walking the process handle while it runs.

Writes a markdown report (default `build/rule_table_footprint.md`) and the
same rows as TSV next to it, so the numbers can be re-analysed without
re-running the sweep.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import signal
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
RM_DRIVER = os.environ.get(
    "RM_DRIVER",
    os.path.join(REPO_ROOT, "build", "release", "rm_driver"),
)
DEFAULT_REPORT = os.path.join(REPO_ROOT, "build", "rule_table_footprint.md")
TIMEOUT_FILE = os.path.join(SCRIPT_DIR, "model_timeouts.tsv")

# ru_maxrss is bytes on the BSDs (macOS included) and kilobytes on Linux.
RSS_UNIT = 1 if sys.platform == "darwin" else 1024


# ---------------------------------------------------------------------------
# Model discovery
#
# The three suites keep their run parameters in different places, which is
# pre-existing: `corpus` and `nfsim_basicmodels` carry a sim_params.tsv next
# to their reference ensembles, while `feature_coverage` takes t_end /
# n_steps from the model's own `simulate` action the way
# benchmark_feature_coverage.py does.  Both are read here rather than
# imported, since neither harness exposes them as a function.
# ---------------------------------------------------------------------------

SUITES = ("feature_coverage", "corpus", "nfsim_basicmodels")


def _params_from_tsv(path: str) -> dict[str, tuple[str, str]]:
    params = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            model = (row.get("#model") or row.get("model") or "").strip()
            if not model or model.startswith("#"):
                continue
            params[model] = (row["t_end"], row["n_steps"])
    return params


def _params_from_bngl(path: str) -> tuple[str, str]:
    """t_end / n_steps from a BNGL `simulate` action, defaulting to 100/100."""
    with open(path) as f:
        text = f.read()
    t_end = n_steps = None
    for args in re.findall(r"simulate\w*\s*\(\s*\{(.*?)\}", text, re.DOTALL):
        m = re.search(r"t_end\s*=>\s*([0-9.eE+\-]+)", args)
        if m:
            t_end = m.group(1)
        m = re.search(r"n_steps\s*=>\s*([0-9]+)", args)
        if m:
            n_steps = m.group(1)
    return t_end or "100", n_steps or "100"


def discover(suite: str) -> list[tuple[str, str, str, str, str]]:
    """[(suite, model, xml_path, t_end, n_steps)] for one suite."""
    out = []
    if suite == "feature_coverage":
        suite_dir = os.path.join(REPO_ROOT, "tests", "models", "feature_coverage")
        xml_dir = os.path.join(suite_dir, "xml")
        for fn in sorted(os.listdir(xml_dir)):
            if not fn.endswith(".xml"):
                continue
            model = fn[:-4]
            bngl = os.path.join(suite_dir, f"{model}.bngl")
            t_end, n_steps = _params_from_bngl(bngl) if os.path.exists(bngl) else ("100", "100")
            out.append((suite, model, os.path.join(xml_dir, fn), t_end, n_steps))
        return out

    ref = "nfsim" if suite == "corpus" else "basicmodels"
    ref_dir = os.path.join(REPO_ROOT, "tests", "reference", ref)
    params = _params_from_tsv(os.path.join(ref_dir, "sim_params.tsv"))
    xml_dir = os.path.join(ref_dir, "xml")
    for model, (t_end, n_steps) in sorted(params.items()):
        xml = os.path.join(xml_dir, f"{model}.xml")
        if os.path.exists(xml):
            out.append((suite, model, xml, t_end, n_steps))
    return out


def load_timeouts() -> dict[str, int]:
    timeouts = {}
    if not os.path.exists(TIMEOUT_FILE):
        return timeouts
    with open(TIMEOUT_FILE) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            model, seconds = line.split("\t")[:2]
            timeouts[model.strip()] = int(seconds)
    return timeouts


# ---------------------------------------------------------------------------
# One run
# ---------------------------------------------------------------------------

TABLES_RE = re.compile(
    r"\[RM tables\] arena=(\d+) rules=(\d+) tabled=(\d+) bytes=(\d+) "
    r"\([\d.]+ MB\) mol=(\d+) sampler=(\d+) reach_bytes=(\d+)"
)
PER_RULE_RE = re.compile(
    r"^  (\S+) \(.*\): rows=(\d+) reach=(\d+) mol=(\d+) sampler=(\d+)$", re.MULTILINE
)


def run_one(xml: str, t_end: str, n_steps: str, timeout: int) -> dict | None:
    """Run rm_driver once and return its footprint row, or None if it failed.

    stdout goes to /dev/null (the .gdat is not what this measures) and
    stderr to a temp file, so nothing has to be drained off a pipe while
    the child runs -- `os.wait4` needs to reap the child itself to get its
    rusage, which `Popen.communicate` would otherwise do first.
    """
    env = {**os.environ, "RM_PRINT_TABLES": "1"}
    with tempfile.TemporaryFile(mode="w+") as errf:
        t0 = time.monotonic()
        proc = subprocess.Popen(
            [RM_DRIVER, xml, t_end, n_steps, "1"],
            stdout=subprocess.DEVNULL,
            stderr=errf,
            env=env,
        )
        deadline = t0 + timeout
        while True:
            pid, status, rusage = os.wait4(proc.pid, os.WNOHANG)
            if pid != 0:
                break
            if time.monotonic() > deadline:
                proc.kill()
                os.wait4(proc.pid, 0)
                return None
            time.sleep(0.02)
        wall_s = time.monotonic() - t0
        proc.returncode = os.waitstatus_to_exitcode(status)
        errf.seek(0)
        stderr_text = errf.read()

    if proc.returncode != 0:
        return None
    m = TABLES_RE.search(stderr_text)
    if not m:
        return None
    arena, rules, tabled, table_bytes, mol_bytes, sampler_bytes, reach_bytes = (
        int(g) for g in m.groups()
    )
    per_rule = PER_RULE_RE.findall(stderr_text)
    widest = max((int(r[3]) // max(int(r[1]), 1) for r in per_rule), default=0)
    return {
        "arena": arena,
        "rules": rules,
        "tabled": tabled,
        "table_bytes": table_bytes,
        "mol_bytes": mol_bytes,
        "sampler_bytes": sampler_bytes,
        "reach_bytes": reach_bytes,
        "widest_row": widest,
        "peak_rss": rusage.ru_maxrss * RSS_UNIT,
        "wall_s": wall_s,
    }


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

COLUMNS = [
    "suite",
    "model",
    "rules",
    "tabled",
    "arena",
    "table_bytes",
    "mol_bytes",
    "sampler_bytes",
    "reach_bytes",
    "widest_row",
    "peak_rss",
    "share",
    "wall_s",
]


def mb(n: float) -> str:
    return f"{n / (1024.0 * 1024.0):.1f}"


def write_report(rows: list[dict], failures: list[tuple[str, str]], path: str) -> None:
    rows = sorted(rows, key=lambda r: r["share"], reverse=True)
    tsv_path = os.path.splitext(path)[0] + ".tsv"
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(tsv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    shares = [r["share"] for r in rows]
    shares_sorted = sorted(shares)
    median = shares_sorted[len(shares_sorted) // 2] if shares_sorted else 0.0
    over_10 = [r for r in rows if r["share"] >= 0.10]
    over_1 = [r for r in rows if r["share"] >= 0.01]
    total_bytes = sum(r["table_bytes"] for r in rows)
    total_reach = sum(r["reach_bytes"] for r in rows)

    with open(path, "w") as f:
        f.write("# Per-rule side-table residency across the corpora (issue #71)\n\n")
        f.write(
            f"{len(rows)} models measured, {len(failures)} skipped.  "
            f"`share` is the fraction of the run's peak RSS held in per-rule "
            f"per-molecule tables and Fenwick samplers.\n\n"
        )
        f.write(
            f"- median share **{median:.2%}**, max **{max(shares, default=0):.2%}**\n"
            f"- models at or above 10% of peak: **{len(over_10)}**\n"
            f"- models at or above 1% of peak: **{len(over_1)}**\n"
            f"- largest arena: **{max((r['arena'] for r in rows), default=0):,}** molecules\n"
            f"- largest table footprint: **{mb(max((r['table_bytes'] for r in rows), default=0))} MB**\n"
            f"- sizing every table by what its rule can reach instead of by the arena "
            f"would hold **{total_reach / total_bytes:.1%}** of the same bytes summed over "
            f"the sweep\n\n"
        )
        f.write(
            "| suite | model | rules | tabled | arena | tables (MB) | mol | sampler "
            "| reach (MB) | peak RSS (MB) | share | wall (s) |\n"
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n"
        )
        for r in rows:
            f.write(
                f"| {r['suite']} | {r['model']} | {r['rules']} | {r['tabled']} "
                f"| {r['arena']:,} | {mb(r['table_bytes'])} | {mb(r['mol_bytes'])} "
                f"| {mb(r['sampler_bytes'])} | {mb(r['reach_bytes'])} | {mb(r['peak_rss'])} "
                f"| {r['share']:.2%} | {r['wall_s']:.1f} |\n"
            )
        if failures:
            f.write("\n## Skipped\n\n")
            for model, why in failures:
                f.write(f"- `{model}`: {why}\n")
    print(f"\nWrote {path}\n      {tsv_path}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("models", nargs="*", help="Model ids to measure; omit for every model")
    p.add_argument("--suite", choices=SUITES, action="append", help="Restrict to a suite")
    p.add_argument("--jobs", type=int, default=1, help="Concurrent runs (default 1)")
    p.add_argument("--output", default=DEFAULT_REPORT, help=f"Report path ({DEFAULT_REPORT})")
    args = p.parse_args()

    if not os.path.exists(RM_DRIVER):
        sys.exit(f"ERROR: {RM_DRIVER} not found. Build it: cmake --build build/release")

    work = []
    for suite in args.suite or SUITES:
        work.extend(discover(suite))
    if args.models:
        wanted = set(args.models)
        work = [w for w in work if w[1] in wanted]
        missing = wanted - {w[1] for w in work}
        if missing:
            sys.exit(f"ERROR: no such model(s): {', '.join(sorted(missing))}")

    timeouts = load_timeouts()
    rows: list[dict] = []
    failures: list[tuple[str, str]] = []
    done = 0

    def measure(item):
        suite, model, xml, t_end, n_steps = item
        return item, run_one(xml, t_end, n_steps, timeouts.get(model, 300))

    def record(item, res):
        nonlocal done
        done += 1
        suite, model, _xml, _t_end, _n_steps = item
        if res is None:
            failures.append((model, "driver failed or timed out"))
            print(f"[{done}/{len(work)}] {suite}/{model}: FAILED", flush=True)
            return
        res.update(suite=suite, model=model)
        res["share"] = res["table_bytes"] / res["peak_rss"] if res["peak_rss"] else 0.0
        rows.append(res)
        print(
            f"[{done}/{len(work)}] {suite}/{model}: {res['rules']} rules, "
            f"arena {res['arena']:,}, tables {mb(res['table_bytes'])} MB of "
            f"{mb(res['peak_rss'])} MB peak ({res['share']:.1%}), {res['wall_s']:.1f}s",
            flush=True,
        )

    t0 = time.monotonic()
    if args.jobs > 1:
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(measure, item): item for item in work}
            for fut in as_completed(futures):
                item, res = fut.result()
                record(item, res)
    else:
        for item in work:
            item, res = measure(item)
            record(item, res)

    print(f"\n{len(rows)} measured, {len(failures)} skipped in {time.monotonic() - t0:.0f}s")
    write_report(rows, failures, args.output)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(128 + signal.SIGINT)
