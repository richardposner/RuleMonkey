# NFsim Reference Data — Provenance

**Last full regeneration:** 2026-04-10
**Generator binary:** `build/NFsim` (local build)
**RM repo commit at regen:** `2de71bf` (`Revert "Fixes to issues #48 and #49"`)
**Reps per model:** 100
**Models:** all 71 in `sim_params.tsv`
**Total regen wall time:** ~4.5 h (6 parallel workers, LPT-scheduled)

## Flag policy

Every reference run invokes NFsim with `-bscb` (which implies `-cb`). This
enforces strict BNGL `+` semantics on BOTH sides of a reversible rule:

- **LHS `+`** — reactants must come from different complexes (no
  intra-complex bimolecular binding).
- **RHS `+`** — products must end up in different complexes (product
  molecularity check, applied only when the rule's product pattern actually
  uses `+`).

Ring closure must be written as an explicit unimolecular rule with `.` on the
LHS (e.g. `A(x).B(y) -> A(x!1).B(y!1)`); ring-opening is the reverse and
uses `.` on the product side, so the product-molecularity check must not
reject it.

RuleMonkey enforces the same semantics by default
(`block_same_complex_binding=true`), so RM and this reference are directly
comparable without flag gymnastics.

Per-model flags from `sim_params.tsv` (e.g. `-utl 5` for rings models) are
applied on top of `-bscb`; `dev/regenerate_nfsim_reference.py` deduplicates
duplicate `-bscb` / `-gml` tokens automatically.

## Generator binary provenance

`build/NFsim` is built from the `src/NFsim_src` tree vendored in this
repository. Relevant fixes present in the binary used for the 2026-04-10
regen:

- **UTL fix** (commit `40e6b93`, upstream `RuleWorld/nfsim#54`) —
  auto-calculated UTL = `max_pattern_size + 1` instead of the buggy
  BNG-2.9.3 value of `max_pattern_size`. This affects every model that
  uses multi-molecule observables or multi-molecule reactant patterns
  without an explicit `-utl` in `sim_params.tsv`.

- **Upstream product-molecularity over-rejection — REVERTED LOCALLY**
  (reverts `72b089d`, `2de71bf`). RuleWorld/nfsim commits `f34dd20` and
  `ca8ad99` added a `checkMolecularity()` path that rejects *every*
  unimolecular unbinding where the molecules remain connected through
  another bond, without checking whether the rule's product pattern uses
  `+` or `.`. This over-rejects valid ring-opening rules (e.g. `A_plus_B_rings`).
  The bug is still present in RuleWorld/nfsim master/dev; we carry the
  revert locally.

**Do NOT regenerate reference data with the BNG-2.9.3 bundled NFsim.**
It lacks the UTL fix and will silently produce wrong values for
multi-molecule patterns.

## Stale reference risk

If any of the following change, the reference data is stale and must be
regenerated:

1. The local `src/NFsim_src` tree changes in a way that affects simulator
   output (rule firing, bond-handling, propensity calculation, output
   formatting).
2. `sim_params.tsv` changes (new model, changed `t_end` / `n_steps`, new
   per-model flags).
3. A new model is added to the benchmark corpus.
4. We accept an upstream NFsim fix/change that alters outputs.

Whenever you regenerate, **update this file** with the new date, commit,
and a one-line note on why the regen was needed.

## What's tracked in git

- `ensemble/*.mean.tsv`, `ensemble/*.std.tsv` — **TRACKED** as of
  2026-04-10. These are the 100-rep aggregates consumed by the benchmark
  and are the whole point of the 4.5h regen. ~8.8 MB across 142 files.
- `ensemble/*.tint.tsv` — **TRACKED** as of 2026-04-10. Per-observable
  time-integral stats (`I_mean`, `I_std`, `n_reps`) computed from the
  replicates, one file per model. Consumed by `benchmark_full.py` to
  compute the `tz_max` verdict metric without needing the full
  replicates/ directory. Produced by `dev/compute_noise_floor.py`.
- `noise_floor.tsv` — **TRACKED** as of 2026-04-10. Per-model
  self-split noise-floor distributions of the primary `max_z` and
  verdict `tz_max` metrics, calibrated at sample sizes n=2/3/10 (matching
  smoke/guard/full tiers). Used by `benchmark_full.py` to derive the
  per-model verdict threshold `T_model = max(5.0, 1.2 × tz_p99)`.
  Produced by `dev/compute_noise_floor.py`; refreshed automatically by
  `regenerate_nfsim_reference.py` at the end of any regen.
- `summary.tsv`, `sim_params.tsv`, `xml/*.xml`, `PROVENANCE.md` —
  tracked.
- `replicates/` — **NOT tracked.** 1.9 GB of raw per-rep trajectories.
  The ensembles are the canonical artifact; replicates are disposable
  scratch kept locally to let us re-aggregate a subset if needed. If
  you need to regen ensembles from replicates, use the per-model
  aggregation logic in `dev/regenerate_nfsim_reference.py`. The
  calibration script (`dev/compute_noise_floor.py`) also needs this
  directory to recompute `noise_floor.tsv` / `*.tint.tsv`, so keep at
  least one full copy locally.

## Regen tooling

- Script: `dev/regenerate_nfsim_reference.py`. At the end of a run it
  automatically invokes `dev/compute_noise_floor.py` on any models it
  regenerated, so `noise_floor.tsv` and the per-model `ensemble/*.tint.tsv`
  stats stay in lock-step with the replicate data. Pass `--skip-noise-floor`
  to bypass the refresh (e.g. if you want to regenerate multiple subsets
  and calibrate once at the end).
- Parallelization: use LPT partitioning with 6 workers; the slow models
  (`pltr`, `mlnr`, `testcase2a`, `rm_tlbr`) should each get their own
  worker to avoid head-of-line blocking.
- Output layout:
  - `ensemble/{model}.mean.tsv`, `ensemble/{model}.std.tsv` — 100-rep
    aggregates (primary artifact consumed by the benchmark).
  - `ensemble/{model}.tint.tsv` — per-observable time-integral stats
    (`I_mean`, `I_std`, `n_reps`) for the `tz_max` verdict metric.
    Produced by `dev/compute_noise_floor.py`.
  - `noise_floor.tsv` — per-model self-split distributions of `max_z`
    and `tz_max` at n=2/3/10, used to derive the per-model verdict
    threshold. Produced by `dev/compute_noise_floor.py`.
  - `replicates/{model}/rep_NNN_{model}.gdat` — raw per-rep trajectories.
  - `summary.tsv` — per-model wall-time summary, used by
    `benchmark_full.py` for progress sorting.

## History

| Date | Commit | Scope | Notes |
|---|---|---|---|
| 2026-04-10 | `2de71bf` | All 71 models | First full regen with `-bscb` + locally-reverted upstream #61 bug |
