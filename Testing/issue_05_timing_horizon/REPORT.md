# Issue #5 — RM vs NFsim time-horizon scaling

**Question (per @wshlavacek, citing RGP):** original RuleMonkey's
per-event cost reportedly grew as the simulation horizon `t_end`
elongated.  Does the new RuleMonkey 3.x exhibit the same artifact?

**Answer: No.** Across 3 mass-action / signaling models and 2
structural-aggregation models, RM's wall-time scales linearly (or
sub-linearly, due to fixed-cost amortization) with event count, and
event count scales linearly with `t_end`.  Per-event time is *lowest*
at the longest horizons we tested — the opposite of the reported
artifact.

NFsim shows the same scaling pattern, so the experiment is also a
sanity check that the test methodology is sound.

---

## Methodology

- **Models** (all from
  [BNGL-Models](https://github.com/wshlavacek/BNGL-Models),
  NFsim-compatible):
  - `Lipniacki2006` — NFkB-like, small mass-action network
  - `Samoilov2005_FutileCycle` — very high event rate (≥ 10⁶ events / time unit)
  - `HBF1998_brusselator` — stiff oscillator
  - `tlbr_macken1982` — Trivalent-Ligand Bivalent-Receptor aggregation
  - `blbr_dembo1978` — Symmetric Bivalent-Ligand Bivalent-Receptor crosslinking
- **Excluded** (NFsim refuses): `Goldbeter1996` (function-of-function
  unsupported), `Mueller2006_RepLeaky_n3` (composite functions of
  observables unsupported), `Bergman1989` (negative-propensity guard
  trips once `Obs_I < I_b`).
- **Ladder**: 4-point geometric, spanning 3 decades of `t_end` per model.
- **Per cell**: 3 reps, distinct seeds; median reported.  `-bscb` on both
  engines.  `n_steps = 400` held constant across the ladder so
  observable-recording cost stays fixed and the comparison isolates
  per-event engine cost.
- **Metric**: wall-time per accepted SSA event.  A *rising* per-event
  time at longer horizons would be the original-RM signature.
- **Harness**: `harness/benchmark_time_horizon.py`.  Raw per-rep data
  in `results/timings.csv`; auto-generated per-model tables in
  `results/summary.md`.

## Headline numbers (median ns/event, short vs long horizon)

| Model | Engine | ns/evt @ short | ns/evt @ long | Drift (long/short) |
|---|---|---:|---:|---:|
| Lipniacki2006              | rm    |   1,502 |   1,363 | 0.91× |
| Lipniacki2006              | nfsim |     821 |     489 | 0.60× |
| Samoilov2005_FutileCycle   | rm    |   3,604 |   3,481 | 0.97× |
| Samoilov2005_FutileCycle   | nfsim |   1,745 |     865 | 0.50× |
| HBF1998_brusselator        | rm    |   1,592 |   1,310 | 0.82× |
| HBF1998_brusselator        | nfsim |     985 |     537 | 0.55× |
| tlbr_macken1982            | rm    | 167,016 |  24,490 | 0.15× |
| tlbr_macken1982            | nfsim | 123,943 |  16,751 | 0.14× |
| blbr_dembo1978             | rm    |  61,351 |  10,804 | 0.18× |
| blbr_dembo1978             | nfsim |  53,901 |   7,990 | 0.15× |

No model on either engine has a drift > 1.0×.  The original-RM
artifact threshold (~1.5× per decade was the rough ballpark RGP
described) would manifest as drift ≫ 1; we see the opposite.

## Why is per-event time *lower* at long horizons?

Two effects, both expected:

1. **Fixed startup cost amortizes**.  Each rep includes process spawn,
   XML parse, RNG seeding, and exit (~30-60 ms on this machine).  At
   `t_end = 20` and 40 K events, that's 1 µs/event of pure overhead;
   at `t_end = 20000` and 40 M events, it's effectively free.
   This is the dominant effect for the mass-action models (drift
   0.5-0.97×).

2. **Aggregation models converge to a steady-state distribution**.
   For TLBR and BLBR, early in the trajectory the network-free engine
   spends time creating and labeling new species as aggregates form.
   Once the cluster-size distribution stabilizes, each event is
   cheaper — same species table, no new labels, hot caches.  This is
   why the aggregation models show the steepest "improvement"
   (drift 0.15-0.18×).  It's the inverse of the artifact under
   investigation.

## RM vs NFsim wall-time

For the record (not the question issue #5 was asking):

| Model | RM wall @ longest horizon | NFsim wall @ longest horizon | RM speedup |
|---|---:|---:|---:|
| Lipniacki2006            | 54.7 s | 19.6 s | 0.36× |
| Samoilov2005_FutileCycle | 76.1 s | 18.9 s | 0.25× |
| HBF1998_brusselator      | 46.2 s | 18.9 s | 0.41× |
| tlbr_macken1982          |  6.3 s |  4.3 s | 0.68× |
| blbr_dembo1978           | 29.4 s | 21.7 s | 0.74× |

RM is currently 1.4× — 4× slower than NFsim wall-time on these models.
The gap is smallest on the aggregation models (RM 0.68-0.74× NFsim),
largest on the high-event-rate Samoilov model (RM 0.25× NFsim).
**This is consistent with the rest of the repo's benchmark history and
is a separate optimization story from the horizon-scaling question.**

## Verdict on issue #5

**RM 3.x scales linearly with the simulation horizon, on both mass-action
and aggregation workloads, across 3 decades of `t_end`.**  Whatever
slowed the original RM as horizons lengthened, it is not present in
the cleanroom rewrite.

## Reproduce

```bash
cd /Users/wish/Code/RuleMonkey
python3 harness/benchmark_time_horizon.py
# Writes Testing/issue_05_timing_horizon/results/{timings.csv,summary.md}
```

Inputs in `Testing/issue_05_timing_horizon/models/` (BNGL),
`Testing/issue_05_timing_horizon/xml/` (NFsim XML, regenerable from
the BNGL via `BNG2.pl writeXML()`).
