#!/usr/bin/env python3
"""RM vs NFsim time-horizon scaling — tests whether RM's per-event cost
grows with t_end.

Context (issue #5): the original RM (NFsim's ancestor) reportedly
slowed down disproportionately as simulation horizons lengthened.  This
script sweeps a geometric ladder of t_end values on a small set of
real-world models, runs N reps per (engine, model, t_end), and reports
wall-time-per-event so a horizon-dependent slowdown would show up as
rising per-event time.

Inputs are fixed in MODELS below: five real-world models from
BNGL-Models that both engines accept — three mass-action / signaling
(Lipniacki2006, Samoilov2005_FutileCycle, HBF1998_brusselator) plus
two structural-aggregation models (TLBR Macken1982, BLBR Dembo1978).
Aggregation models stress the network-free path, which historically
slowed disproportionately at long horizons in the original RM.  Each
has a hand-picked geometric ladder spanning ~3 decades of t_end with
the longest rep capped at ~120s wall.

Models that NFsim refuses are excluded: Goldbeter1996 and
Mueller2006_RepLeaky_n3 use BNGL function features NFsim doesn't
implement; Bergman1989 trips NFsim's negative-propensity guard once
its `Obs_I` observable drops below the `I_b` setpoint.  This is a
comparison benchmark and RM-only rows don't serve the question.

n_steps is held constant across the ladder so observable-recording cost
stays fixed; the comparison isolates per-event engine cost from output
cost.

Outputs:
  Testing/issue_05_timing_horizon/results/timings.csv      (per-rep)
  Testing/issue_05_timing_horizon/results/summary.md       (writeup)
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
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
TEST_DIR = REPO_ROOT / "Testing" / "issue_05_timing_horizon"
XML_DIR = TEST_DIR / "xml"
RESULTS_DIR = TEST_DIR / "results"

RM_DRIVER = Path(os.environ.get("RM_DRIVER", REPO_ROOT / "build" / "release" / "rm_driver"))
NFSIM_BIN = Path(os.environ.get("NFSIM_BIN", "/Users/wish/Simulations/BioNetGen-2.9.3/bin/NFsim"))

# (model_stem, ladder of t_end values).  Ladders chosen per-model to keep
# the longest rep ~30-120s wall while spanning >= 3 decades of t_end.
MODELS: list[tuple[str, list[float]]] = [
    ("Lipniacki2006", [20.0, 200.0, 2000.0, 20000.0]),
    ("Samoilov2005_FutileCycle", [0.018, 0.18, 1.8, 18.0]),
    ("HBF1998_brusselator", [10.0, 100.0, 1000.0, 10000.0]),
    ("tlbr_macken1982", [30.0, 300.0, 3000.0, 30000.0]),
    ("blbr_dembo1978", [30.0, 300.0, 3000.0, 30000.0]),
]

N_STEPS = 400  # constant across ladder to isolate per-event cost
N_REPS = 3
TIMEOUT_S = 300

RM_TIMING_RE = re.compile(r"events=(\d+)\s+null=(\d+)\s+total=([\d.]+)s\s+wall=([\d.]+)s")
NFSIM_EVENTS_RE = re.compile(r"You just simulated\s+(\d+)\s+reactions")


@dataclass
class Rep:
    engine: str
    model: str
    t_end: float
    seed: int
    wall_s: float
    events: int
    error: str = ""


def run_rm(xml: Path, t_end: float, seed: int) -> Rep:
    env = dict(os.environ)
    env["RM_PRINT_TIMING"] = "1"
    cmd = [str(RM_DRIVER), str(xml), str(t_end), str(N_STEPS), str(seed)]
    t0 = time.monotonic()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT_S, env=env)
    except subprocess.TimeoutExpired:
        return Rep("rm", xml.stem, t_end, seed, -1.0, -1, f"timeout>{TIMEOUT_S}s")
    wall = time.monotonic() - t0
    if r.returncode != 0:
        return Rep("rm", xml.stem, t_end, seed, wall, -1, f"rc={r.returncode}: {r.stderr[:120]}")
    m = RM_TIMING_RE.search(r.stderr)
    events = int(m.group(1)) if m else -1
    return Rep("rm", xml.stem, t_end, seed, wall, events)


def run_nfsim(xml: Path, t_end: float, seed: int) -> Rep:
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "out.gdat"
        cmd = [
            str(NFSIM_BIN),
            "-xml",
            str(xml),
            "-sim",
            str(t_end),
            "-oSteps",
            str(N_STEPS),
            "-seed",
            str(seed),
            "-bscb",
            "-o",
            str(out),
        ]
        t0 = time.monotonic()
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT_S, cwd=tmp)
        except subprocess.TimeoutExpired:
            return Rep("nfsim", xml.stem, t_end, seed, -1.0, -1, f"timeout>{TIMEOUT_S}s")
        wall = time.monotonic() - t0
        if r.returncode != 0 or not out.exists() or out.stat().st_size == 0:
            return Rep("nfsim", xml.stem, t_end, seed, wall, -1, f"rc={r.returncode}")
        m = NFSIM_EVENTS_RE.search(r.stdout)
        events = int(m.group(1)) if m else -1
        return Rep("nfsim", xml.stem, t_end, seed, wall, events)


def fmt_ns_per_event(wall_s: float, events: int) -> str:
    if events <= 0 or wall_s < 0:
        return "—"
    ns = (wall_s / events) * 1e9
    if ns < 1000:
        return f"{ns:.0f} ns"
    if ns < 1e6:
        return f"{ns / 1000:.1f} µs"
    return f"{ns / 1e6:.2f} ms"


def fmt_wall(wall_s: float) -> str:
    if wall_s < 0:
        return "—"
    if wall_s < 0.1:
        return f"{wall_s * 1000:.1f} ms"
    if wall_s < 10:
        return f"{wall_s * 1000:.0f} ms"
    return f"{wall_s:.1f} s"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reps", type=int, default=N_REPS, help=f"default {N_REPS}")
    ap.add_argument(
        "--quick",
        action="store_true",
        help="run only the shortest 2 t_end points per model (smoke test)",
    )
    args = ap.parse_args()

    if not RM_DRIVER.exists():
        sys.exit(f"RM_DRIVER missing: {RM_DRIVER}")
    if not NFSIM_BIN.exists():
        sys.exit(f"NFSIM_BIN missing: {NFSIM_BIN}")
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    reps: list[Rep] = []
    t_overall = time.monotonic()
    for model, ladder in MODELS:
        xml = XML_DIR / f"{model}.xml"
        if not xml.exists():
            print(f"WARN: missing XML for {model}; skipping", file=sys.stderr)
            continue
        sweep = ladder[:2] if args.quick else ladder
        for t_end in sweep:
            for engine, runner in (("rm", run_rm), ("nfsim", run_nfsim)):
                for rep_i in range(args.reps):
                    seed = rep_i + 1
                    rep = runner(xml, t_end, seed)
                    reps.append(rep)
                    elapsed = time.monotonic() - t_overall
                    if rep.error:
                        msg = rep.error
                    else:
                        msg = (
                            f"wall={fmt_wall(rep.wall_s)} events={rep.events} "
                            f"per_evt={fmt_ns_per_event(rep.wall_s, rep.events)}"
                        )
                    print(
                        f"[{elapsed:6.1f}s] {model:<28} t_end={t_end:<9g} "
                        f"{engine:<6} seed={seed}: {msg}",
                        file=sys.stderr,
                    )

    # CSV (per-rep)
    csv_path = RESULTS_DIR / "timings.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["engine", "model", "t_end", "seed", "wall_s", "events", "error"])
        for r in reps:
            w.writerow([r.engine, r.model, r.t_end, r.seed, f"{r.wall_s:.6f}", r.events, r.error])

    # Aggregate per (engine, model, t_end)
    @dataclass
    class Cell:
        wall_med: float
        events_med: float
        n_ok: int
        n_err: int

    def aggregate(rows: list[Rep]) -> Cell | None:
        walls = [r.wall_s for r in rows if r.events > 0 and r.wall_s > 0]
        events = [r.events for r in rows if r.events > 0]
        if not walls:
            return None
        return Cell(
            wall_med=statistics.median(walls),
            events_med=statistics.median(events),
            n_ok=len(walls),
            n_err=len(rows) - len(walls),
        )

    agg: dict[tuple[str, str, float], Cell] = {}
    grouped: dict[tuple[str, str, float], list[Rep]] = {}
    for r in reps:
        grouped.setdefault((r.engine, r.model, r.t_end), []).append(r)
    for key, rows in grouped.items():
        c = aggregate(rows)
        if c is not None:
            agg[key] = c

    # Write summary markdown
    md_path = RESULTS_DIR / "summary.md"
    with md_path.open("w") as f:
        f.write("# RM vs NFsim — time-horizon scaling (issue #5)\n\n")
        f.write(
            "Tests whether RuleMonkey 3.x exhibits horizon-dependent slowdown "
            "(the artifact RGP observed in original RM).  Each row is the median "
            f"of {args.reps} reps; `n_steps={N_STEPS}` held constant across the "
            "ladder so observable-recording cost stays fixed and the comparison "
            "isolates per-event engine cost.\n\n"
        )
        f.write(
            "**Reading the table:** if per-event time stays roughly constant down "
            "each model's column, the engine scales linearly with horizon (no "
            "degradation).  A *rising* per-event time at longer horizons would be "
            "the original-RM signature.\n\n"
        )
        f.write(
            f"Generated: `{time.strftime('%Y-%m-%d %H:%M:%S %Z')}` "
            f"(RuleMonkey build: `{RM_DRIVER}`)\n\n"
        )

        for model, ladder in MODELS:
            f.write(f"## `{model}`\n\n")
            f.write(
                "| t_end | RM wall | RM events | RM ns/event | "
                "NFsim wall | NFsim events | NFsim ns/event | "
                "RM speedup |\n"
            )
            f.write("|---:|---:|---:|---:|---:|---:|---:|---:|\n")
            ladder_used = ladder[:2] if args.quick else ladder
            for t_end in ladder_used:
                rm = agg.get(("rm", model, t_end))
                nf = agg.get(("nfsim", model, t_end))
                rm_wall = fmt_wall(rm.wall_med) if rm else "—"
                rm_evt = f"{int(rm.events_med):,}" if rm else "—"
                rm_ns = fmt_ns_per_event(rm.wall_med, int(rm.events_med)) if rm else "—"
                nf_wall = fmt_wall(nf.wall_med) if nf else "—"
                nf_evt = f"{int(nf.events_med):,}" if nf else "—"
                nf_ns = fmt_ns_per_event(nf.wall_med, int(nf.events_med)) if nf else "—"
                if rm and nf and rm.wall_med > 0:
                    sp = f"{nf.wall_med / rm.wall_med:.2f}×"
                else:
                    sp = "—"
                f.write(
                    f"| {t_end:g} | {rm_wall} | {rm_evt} | {rm_ns} | "
                    f"{nf_wall} | {nf_evt} | {nf_ns} | {sp} |\n"
                )
            f.write("\n")

        # Horizon-scaling check: for each engine+model, compare ns/event at
        # smallest vs largest t_end.  Ratio > 1.5× flags potential degradation.
        f.write("## Per-event time: smallest vs largest horizon\n\n")
        f.write(
            "Drift > 1.5× (slower at longer horizon) would suggest engine cost "
            "growing with simulation state or time, not just with event count.\n\n"
        )
        f.write("| Model | Engine | ns/evt @ short | ns/evt @ long | Drift (long/short) |\n")
        f.write("|---|---|---:|---:|---:|\n")
        for model, ladder in MODELS:
            ladder_used = ladder[:2] if args.quick else ladder
            short_t, long_t = ladder_used[0], ladder_used[-1]
            for engine in ("rm", "nfsim"):
                a = agg.get((engine, model, short_t))
                b = agg.get((engine, model, long_t))
                if not a or not b or a.events_med <= 0 or b.events_med <= 0:
                    continue
                ns_a = (a.wall_med / a.events_med) * 1e9
                ns_b = (b.wall_med / b.events_med) * 1e9
                drift = ns_b / ns_a if ns_a > 0 else float("nan")
                flag = " ⚠" if drift > 1.5 else ""
                f.write(
                    f"| `{model}` | {engine} | "
                    f"{ns_a:.0f} ns | {ns_b:.0f} ns | "
                    f"{drift:.2f}×{flag} |\n"
                )
        f.write("\n")

        # Quick interpretation footer
        f.write("## Notes on methodology\n\n")
        f.write(
            "- Wall time is end-to-end subprocess wall (process spawn + XML load + "
            "simulate + observable record + exit).  Same for both engines, so "
            "ratios are fair; absolute numbers shift per machine.\n"
            "- RM event count comes from `RM_PRINT_TIMING=1` stderr "
            "(`events=N null=N`); NFsim event count from stdout "
            "(`You just simulated N reactions`).  Both count accepted SSA events.\n"
            "- `-bscb` (Bind-Single-Copy-Bond) enabled on both engines.\n"
            "- `n_steps` held constant; per-event cost should track engine work, "
            "not observable I/O.\n"
            "- Each rep uses a distinct seed (1..N), so event totals differ "
            "between reps; medians are reported.\n"
        )

    print(f"\nWrote {csv_path}\nWrote {md_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
