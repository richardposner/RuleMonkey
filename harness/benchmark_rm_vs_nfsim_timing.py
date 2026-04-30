#!/usr/bin/env python3
"""RuleMonkey vs NFsim wall-time comparison.

Runs every model from all three suites (corpus, feature_coverage,
nfsim_basicmodels) sequentially in both engines, captures wall time
across N reps each, and writes a sorted comparison table to
``docs/timing_comparison.md``.

Methodology:
  * 3 reps per engine per model.  Median reported (robust to noise on
    small N); spread shown via min/max.
  * Sequential runs only — no parallel workers.  Parallel runs would
    introduce scheduler noise that obscures the comparison.
  * Same XML, t_end, n_steps for both engines.  Seeds 1..3.
  * Wall time is end-to-end subprocess time including process startup
    and XML load — measured identically for both engines, so the
    comparison is fair.
  * NFsim binary path read from ``NFSIM_BIN`` env var (default
    ``~/Code/nfsim-rm/build/NFsim``); RM driver from ``RM_DRIVER``.
  * Models where NFsim refuses to load or produces no output are
    marked ``N/A`` with the captured stderr summary as the reason.
    NFsim that runs but produces wrong observables (e.g.
    ``ft_tfun`` — TFUN handler returns zero rate) is still included
    in the timing table; the artifact is about efficiency, not
    correctness.

Output: ``docs/timing_comparison.md`` (committed).
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import statistics
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
RM_DRIVER = Path(os.environ.get("RM_DRIVER", REPO_ROOT / "build" / "release" / "rm_driver"))
NFSIM_BIN = Path(os.environ.get("NFSIM_BIN", Path.home() / "Code" / "nfsim-rm" / "build" / "NFsim"))

DEFAULT_N_REPS = 3
DEFAULT_TIMEOUT_S = 300  # per single rep

OUT_PATH = REPO_ROOT / "docs" / "timing_comparison.md"


@dataclass
class Suite:
    name: str
    xml_dir: Path
    sim_params_tsv: Path
    bngl_dir: Path | None = None  # for fallback parsing when TSV is incomplete
    flags_col: str = "nfsim_flags"


SUITES: list[Suite] = [
    Suite(
        name="corpus",
        xml_dir=REPO_ROOT / "tests" / "reference" / "nfsim" / "xml",
        sim_params_tsv=REPO_ROOT / "tests" / "reference" / "nfsim" / "sim_params.tsv",
    ),
    Suite(
        name="feature_coverage",
        xml_dir=REPO_ROOT / "tests" / "models" / "feature_coverage" / "xml",
        sim_params_tsv=REPO_ROOT / "tests" / "models" / "feature_coverage" / "sim_params.tsv",
        bngl_dir=REPO_ROOT / "tests" / "models" / "feature_coverage",
    ),
    Suite(
        name="nfsim_basicmodels",
        xml_dir=REPO_ROOT / "tests" / "reference" / "basicmodels" / "xml",
        sim_params_tsv=REPO_ROOT / "tests" / "reference" / "basicmodels" / "sim_params.tsv",
    ),
]


def _parse_simulate_action(bngl_path: Path) -> tuple[float, int, list[str]] | None:
    """Read t_end, n_steps, NFsim flags from a BNGL's simulate() action.

    Mirrors benchmark_feature_coverage.extract_sim_params.
    Returns None when no simulate action is found.
    """
    if not bngl_path.exists():
        return None
    text = bngl_path.read_text()
    t_end = None
    n_steps = None
    flags: list[str] = []
    for m in re.finditer(r"simulate\w*\s*\(\s*\{([^}]+)\}", text, flags=re.DOTALL):
        args = m.group(1)
        tm = re.search(r"t_end\s*=>\s*([0-9.eE+\-]+)", args)
        nm = re.search(r"n_steps\s*=>\s*([0-9]+)", args)
        if tm:
            t_end = float(tm.group(1))
        if nm:
            n_steps = int(nm.group(1))
        pm = re.search(r'param\s*=>\s*"([^"]*)"', args)
        if pm:
            for tok in pm.group(1).split():
                flags.append(tok)
    if t_end is None or n_steps is None:
        return None
    return t_end, n_steps, flags


@dataclass
class ModelEntry:
    suite: str
    model: str
    xml: Path
    t_end: float
    n_steps: int
    nfsim_flags: list[str] = field(default_factory=list)


def discover_models(suite: Suite) -> list[ModelEntry]:
    """Locate every XML in the suite, with sim params from TSV or (fallback) BNGL.

    Some legacy TSVs use ``#model`` as the model column name; that's
    accepted alongside ``model``.  The BNGL fallback is needed for
    feature_coverage where the TSV lags behind the model directory.
    """
    # Build a map model -> (t_end, n_steps, flags) from the TSV.
    tsv_params: dict[str, tuple[float, int, list[str]]] = {}
    if suite.sim_params_tsv.exists():
        with suite.sim_params_tsv.open() as f:
            reader = csv.DictReader(f, delimiter="\t")
            # Some TSVs have a header row beginning with `#model`; csv.DictReader
            # keeps the leading `#` in the column name.
            model_col = "model" if "model" in (reader.fieldnames or []) else "#model"
            for row in reader:
                m = (row.get(model_col) or "").strip()
                if not m:
                    continue
                try:
                    t_end = float(row.get("t_end", "0"))
                    n_steps = int(row.get("n_steps", "0"))
                except ValueError:
                    continue
                flags_raw = row.get(suite.flags_col, "") or ""
                tsv_params[m] = (t_end, n_steps, flags_raw.split())
    else:
        print(f"WARN: missing sim_params for {suite.name}: {suite.sim_params_tsv}", file=sys.stderr)

    # Enumerate XMLs, look up params with BNGL fallback.
    out: list[ModelEntry] = []
    if not suite.xml_dir.exists():
        return out
    for xml in sorted(suite.xml_dir.glob("*.xml")):
        model = xml.stem
        if model in tsv_params:
            t_end, n_steps, flags = tsv_params[model]
        elif suite.bngl_dir is not None:
            parsed = _parse_simulate_action(suite.bngl_dir / f"{model}.bngl")
            if parsed is None:
                continue
            t_end, n_steps, flags = parsed
        else:
            continue
        out.append(ModelEntry(suite.name, model, xml, t_end, n_steps, flags))
    return out


@dataclass
class TimingResult:
    walls: list[float] = field(default_factory=list)  # seconds, per rep
    error: str = ""  # populated on engine failure

    @property
    def ok(self) -> bool:
        return bool(self.walls) and not self.error

    def median(self) -> float:
        return statistics.median(self.walls) if self.walls else float("nan")

    def min(self) -> float:
        return min(self.walls) if self.walls else float("nan")

    def max(self) -> float:
        return max(self.walls) if self.walls else float("nan")


def _classify_nfsim_error(stderr: str, stdout: str) -> str:
    blob = (stderr + "\n" + stdout).lower()
    if "couldn't create a system" in blob or "i don't know what you did" in blob:
        return "NFsim refused XML"
    if "func factory" in blob:
        return "NFsim FuncFactory error (unsupported expression)"
    if "type 'sat'" in blob:
        return "NFsim refuses Sat type"
    if "tfun" in blob and ("error" in blob or "quitting" in blob):
        return "NFsim TFUN handler error"
    return "NFsim produced no output"


def run_rm(entry: ModelEntry, seed: int, timeout: float) -> tuple[float, str]:
    """Returns (wall_seconds, error_msg). wall < 0 on error."""
    cmd = [str(RM_DRIVER), str(entry.xml), str(entry.t_end), str(entry.n_steps), str(seed)]
    if "-bscb" in entry.nfsim_flags:
        # RM defaults to bscb=true, but make it explicit so the matrix is consistent.
        pass
    elif "-no-bscb" in entry.nfsim_flags:
        cmd.append("-no-bscb")
    t0 = time.monotonic()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return -1.0, f"RM timeout >{timeout}s"
    wall = time.monotonic() - t0
    if r.returncode != 0:
        return -1.0, f"RM rc={r.returncode}: {r.stderr.strip()[:80]}"
    return wall, ""


def run_nfsim(entry: ModelEntry, seed: int, timeout: float) -> tuple[float, str]:
    """Returns (wall_seconds, error_msg). wall < 0 on error."""
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "out.gdat"
        # Stage tfun files alongside in case of file-relative lookups.
        for tfun in entry.xml.parent.glob("*.tfun"):
            (Path(tmp) / tfun.name).write_bytes(tfun.read_bytes())
        # Also try parent (BNGL-side) directory.
        bngl_dir = entry.xml.parent.parent
        for tfun in bngl_dir.glob("*.tfun"):
            tgt = Path(tmp) / tfun.name
            if not tgt.exists():
                tgt.write_bytes(tfun.read_bytes())
        cmd = [
            str(NFSIM_BIN),
            "-xml",
            str(entry.xml),
            "-sim",
            str(entry.t_end),
            "-oSteps",
            str(entry.n_steps),
            "-seed",
            str(seed),
            "-bscb",
            "-o",
            str(out),
        ]
        # Drop nfsim_flags entries that duplicate cmd-line flags we already
        # added (script controls sim time, output, seed, bscb explicitly).
        # nfsim_basicmodels TSV bakes "-sim N -oSteps M -bscb" into its
        # nfsim_flags column; passing it through verbatim produces duplicate
        # flags that NFsim rejects.
        skip_with_arg = {"-sim", "-oSteps", "-seed", "-o", "-xml"}
        skip_solo = {"-bscb"}
        i = 0
        while i < len(entry.nfsim_flags):
            flag = entry.nfsim_flags[i]
            if flag in skip_with_arg:
                i += 2  # also skip the arg
                continue
            if flag in skip_solo:
                i += 1
                continue
            cmd.append(flag)
            i += 1
        t0 = time.monotonic()
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, cwd=tmp)
        except subprocess.TimeoutExpired:
            return -1.0, f"NFsim timeout >{timeout}s"
        wall = time.monotonic() - t0
        if r.returncode != 0 or not out.exists() or out.stat().st_size == 0:
            return -1.0, _classify_nfsim_error(r.stderr, r.stdout)
        return wall, ""


def time_engine(entry: ModelEntry, runner, n_reps: int, timeout: float) -> TimingResult:
    res = TimingResult()
    for rep in range(n_reps):
        seed = rep + 1
        wall, err = runner(entry, seed, timeout)
        if err:
            # Stop on first error and report it; partial timings don't combine
            # cleanly with N/A models.
            res.error = err
            res.walls.clear()
            return res
        res.walls.append(wall)
    return res


def fmt_ms(seconds: float) -> str:
    """Render seconds as ms, with adaptive precision."""
    ms = seconds * 1000.0
    if ms < 1:
        return f"{ms:.2f} ms"
    if ms < 100:
        return f"{ms:.1f} ms"
    if ms < 10_000:
        return f"{ms:.0f} ms"
    return f"{seconds:.1f} s"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--reps",
        type=int,
        default=DEFAULT_N_REPS,
        help=f"reps per engine per model (default {DEFAULT_N_REPS})",
    )
    ap.add_argument(
        "--timeout",
        type=float,
        default=DEFAULT_TIMEOUT_S,
        help=f"per-rep timeout in seconds (default {DEFAULT_TIMEOUT_S})",
    )
    ap.add_argument(
        "--out", type=Path, default=OUT_PATH, help=f"output markdown path (default {OUT_PATH})"
    )
    ap.add_argument(
        "--limit", type=int, default=None, help="run only the first N models (for smoke testing)"
    )
    ap.add_argument("--suite", action="append", help="restrict to named suite(s); default: all")
    args = ap.parse_args()

    # Sanity-check binaries
    if not RM_DRIVER.exists():
        sys.exit(f"ERROR: RM_DRIVER not found at {RM_DRIVER}")
    if not NFSIM_BIN.exists():
        sys.exit(f"ERROR: NFSIM_BIN not found at {NFSIM_BIN}")

    suites = SUITES if not args.suite else [s for s in SUITES if s.name in args.suite]
    entries: list[ModelEntry] = []
    for s in suites:
        suite_entries = discover_models(s)
        print(f"# {s.name}: {len(suite_entries)} models", file=sys.stderr)
        entries.extend(suite_entries)
    if args.limit:
        entries = entries[: args.limit]

    print(f"# Total: {len(entries)} models, {args.reps} reps each, both engines", file=sys.stderr)

    # Run sequentially.
    rows = []
    t_start = time.monotonic()
    for i, e in enumerate(entries, 1):
        rm = time_engine(e, run_rm, args.reps, args.timeout)
        nf = time_engine(e, run_nfsim, args.reps, args.timeout)
        if rm.ok and nf.ok:
            speedup = nf.median() / rm.median() if rm.median() > 0 else float("inf")
            sp_str = f"{speedup:.2f}×"
        elif rm.ok and not nf.ok:
            sp_str = "—"
        elif not rm.ok and nf.ok:
            sp_str = "—"
        else:
            sp_str = "—"
        rows.append(
            {
                "suite": e.suite,
                "model": e.model,
                "t_end": e.t_end,
                "n_steps": e.n_steps,
                "rm": rm,
                "nf": nf,
                "speedup_str": sp_str,
            }
        )
        elapsed = time.monotonic() - t_start
        print(
            f"[{i:>3}/{len(entries)}] {e.suite:<18} {e.model:<40} "
            f"RM={fmt_ms(rm.median()) if rm.ok else 'ERR':>10}  "
            f"NF={fmt_ms(nf.median()) if nf.ok else 'N/A':>10}  "
            f"sp={sp_str:>8}  ({elapsed:.0f}s)",
            file=sys.stderr,
        )

    # Sort by RM median time descending — the "expensive" models lead.
    rows.sort(key=lambda r: r["rm"].median() if r["rm"].ok else -1, reverse=True)

    # Compute summary stats over the both-OK subset.
    speedups = []
    for r in rows:
        if r["rm"].ok and r["nf"].ok and r["rm"].median() > 0:
            speedups.append(r["nf"].median() / r["rm"].median())
    n_total = len(rows)
    n_both_ok = len(speedups)
    n_nf_na = sum(1 for r in rows if r["rm"].ok and not r["nf"].ok)
    n_rm_err = sum(1 for r in rows if not r["rm"].ok)
    median_speedup = statistics.median(speedups) if speedups else float("nan")
    geom_mean = statistics.geometric_mean(speedups) if speedups else float("nan")
    rm_faster = sum(1 for s in speedups if s > 1.0)
    nf_faster = sum(1 for s in speedups if s < 1.0)

    # Write markdown.
    out_path: Path = args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        f.write("# RuleMonkey vs NFsim — wall-time comparison\n\n")
        f.write(
            "Single-machine, sequential, end-to-end subprocess wall time.  "
            f"{args.reps} reps per engine per model; median reported.  "
            "Both engines invoked from a fresh process per rep (no warm caches).  "
            "Same XML, same t_end, same n_steps, same seed sequence.\n\n"
        )
        f.write(
            "**This is an efficiency report, not a correctness report.**  Models "
            "where NFsim runs but produces incorrect observables (e.g. `ft_tfun` "
            "due to NFsim's TFUN handler returning zero rate) are still included "
            "in the timing table; correctness is covered separately by "
            "`feature_coverage`, `benchmark_full`, and `basicmodels` suites.\n\n"
        )
        f.write(
            f"Generated: `{time.strftime('%Y-%m-%d %H:%M:%S %Z')}` "
            f"(`{__file__.split('/')[-1]}`)\n\n"
        )
        f.write("## Summary\n\n")
        f.write(
            f"- Models scored: **{n_both_ok} / {n_total}** "
            f"(both engines completed all {args.reps} reps)\n"
        )
        f.write(f"- NFsim N/A: **{n_nf_na}** (RM ran; NFsim refused / errored)\n")
        if n_rm_err:
            f.write(f"- RM errored: **{n_rm_err}** (investigate before publishing)\n")
        if speedups:
            f.write(f"- **Median speedup (NFsim wall ÷ RM wall): {median_speedup:.2f}×**\n")
            f.write(f"- Geometric mean speedup: {geom_mean:.2f}×\n")
            f.write(
                f"- RM faster than NFsim on **{rm_faster} / {n_both_ok}** models; "
                f"NFsim faster on {nf_faster}\n"
            )
        f.write("\n")
        f.write("## Caveats\n\n")
        f.write(
            "- Wall time, not CPU time.  Includes process startup, XML load, "
            "and result write — same for both engines, so the ratio is fair, "
            "but absolute numbers shift across machines.\n"
        )
        f.write(
            "- Speedup depends on workload: long-horizon, high-event-rate models "
            "show RM's biggest wins; short-horizon trivial models are near parity "
            "or NFsim-favored due to RM's slightly heavier per-process startup.\n"
        )
        f.write(
            "- Single-replicate timings can be noisy.  See min/max columns for "
            "spread; the median column is what to quote.\n"
        )
        f.write(
            "- N/A in the NFsim column means NFsim refused or errored on at "
            "least one of the 3 reps.  The reason column gives the captured "
            "diagnostic.\n\n"
        )
        f.write("## Per-model results\n\n")
        f.write("Sorted by RM median wall time (most expensive first).\n\n")
        f.write(
            "| # | Suite | Model | t_end | n_steps | RM median | RM min/max | "
            "NFsim median | NFsim min/max | Speedup | Notes |\n"
        )
        f.write("|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|\n")
        for i, r in enumerate(rows, 1):
            rm = r["rm"]
            nf = r["nf"]
            rm_med = fmt_ms(rm.median()) if rm.ok else "—"
            rm_range = f"{fmt_ms(rm.min())} / {fmt_ms(rm.max())}" if rm.ok else "—"
            nf_med = fmt_ms(nf.median()) if nf.ok else "N/A"
            nf_range = f"{fmt_ms(nf.min())} / {fmt_ms(nf.max())}" if nf.ok else "—"
            notes = ""
            if not rm.ok:
                notes = f"RM error: {rm.error}"
            elif not nf.ok:
                notes = nf.error
            f.write(
                f"| {i} | {r['suite']} | `{r['model']}` | "
                f"{r['t_end']:g} | {r['n_steps']} | "
                f"{rm_med} | {rm_range} | {nf_med} | {nf_range} | "
                f"{r['speedup_str']} | {notes} |\n"
            )
    print(f"\n# Wrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
