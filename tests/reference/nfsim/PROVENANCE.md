# Reference Data — Provenance

The reference trajectories under `tests/reference/nfsim/` are **mostly** the
output of a gold-standard NFsim build, with a documented set of per-model
exceptions where NFsim was found to produce incorrect output and the
reference was replaced with hand-rolled Gillespie SSA data. See the
"Per-model exceptions" section below.

The directory keeps the `nfsim/` name for historical continuity; consumers
should treat its contents as the canonical correctness reference for the
RuleMonkey benchmark, not as pure NFsim output.

**Last full regeneration:** 2026-04-10
**Generator binary:** `build/NFsim` (local build, in the `nfsim-rm`
development repository — not in this repo)
**Source-repo commit at regen:** `2de71bf` (`Revert "Fixes to issues #48 and #49"`)
**Reps per model:** 100
**Models:** all 71 in `sim_params.tsv`
**Total regen wall time:** ~4.5 h (6 parallel workers, LPT-scheduled)

## Per-model exceptions: SSA-generated reference

The following models were found to produce incorrect output under NFsim,
and their reference data was regenerated with hand-rolled Gillespie SSA
implementations (100 reps each, same `t_end` / `n_steps` as the NFsim
reference). For these models, the ensemble files in `ensemble/` are the
SSA output, **not** NFsim output.

| Model | Replaced | Commit | NFsim error | SSA script |
|---|---|---|---|---|
| `toy_jim` | 2026-04-13 | `9cfe36f` | NFsim over-counts `K(Y~P)` by ~12% on the disjoint transphosphorylation / multi-molecule dephosphorylation patterns. SSA confirmed RM is correct. | [`harness/ssa/toy_jim_ssa.py`](../../../harness/ssa/toy_jim_ssa.py) |
| `rm_tlbr_rings` | 2026-04-12 | `05a9a12` | NFsim under-counts on disjoint + symmetric ring-closure patterns. SSA confirmed RM is correct. | [`harness/ssa/rm_tlbr_rings_ssa.py`](../../../harness/ssa/rm_tlbr_rings_ssa.py) |

The SSA scripts in `harness/ssa/` are sufficient to regenerate these two
models' references from scratch and should be re-run if the corresponding
model definitions in `tests/models/corpus/` change.

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
  replicates/ directory. Produced by the `dev/compute_noise_floor.py`
  script in the `nfsim-rm` sibling repo (see "Regen tooling").
- `noise_floor.tsv` — **TRACKED** as of 2026-04-10. Per-model
  self-split noise-floor distributions of the primary `max_z` and
  verdict `tz_max` metrics, calibrated at sample sizes n=2/3/10 (matching
  smoke/guard/full tiers). Used by `benchmark_full.py` to derive the
  per-model verdict threshold `T_model = max(5.0, 1.2 × tz_p99)`.
  Produced by the `nfsim-rm` repo's `dev/compute_noise_floor.py`;
  refreshed automatically by that repo's `regenerate_nfsim_reference.py`
  at the end of any regen.
- `summary.tsv`, `sim_params.tsv`, `xml/*.xml`, `PROVENANCE.md` —
  tracked.
- `replicates/` — **NOT tracked.** 1.9 GB of raw per-rep trajectories.
  The ensembles are the canonical artifact; replicates are disposable
  scratch kept locally to let us re-aggregate a subset if needed. If
  you need to regen ensembles from replicates, use the per-model
  aggregation logic in the `nfsim-rm` repo's
  `dev/regenerate_nfsim_reference.py`. The calibration script
  (`nfsim-rm`'s `dev/compute_noise_floor.py`) also needs this
  directory to recompute `noise_floor.tsv` / `*.tint.tsv`, so keep at
  least one full copy of the replicates in that repo's tree.

## Regen tooling

- **NFsim reference (69 of 71 models):** regeneration is performed in the
  `nfsim-rm` development repository, which retains the NFsim build
  infrastructure that the cleanroom RuleMonkey repo intentionally does
  not vendor. The driver script is `dev/regenerate_nfsim_reference.py`
  in that repo. At the end of a run it automatically invokes
  `dev/compute_noise_floor.py` on any models it regenerated, so
  `noise_floor.tsv` and the per-model `ensemble/*.tint.tsv` stats stay
  in lock-step with the replicate data. Pass `--skip-noise-floor` to
  bypass the refresh (e.g. if you want to regenerate multiple subsets
  and calibrate once at the end).
- **SSA reference (`toy_jim`, `rm_tlbr_rings`):** regeneration is
  self-contained in this repo; run the corresponding script in
  `harness/ssa/`.
- Parallelization: use LPT partitioning with 6 workers; the slow models
  (`pltr`, `mlnr`, `testcase2a`, `rm_tlbr`) should each get their own
  worker to avoid head-of-line blocking.
- Output layout:
  - `ensemble/{model}.mean.tsv`, `ensemble/{model}.std.tsv` — 100-rep
    aggregates (primary artifact consumed by the benchmark).
  - `ensemble/{model}.tint.tsv` — per-observable time-integral stats
    (`I_mean`, `I_std`, `n_reps`) for the `tz_max` verdict metric.
    Produced by `nfsim-rm`'s `dev/compute_noise_floor.py`.
  - `noise_floor.tsv` — per-model self-split distributions of `max_z`
    and `tz_max` at n=2/3/10, used to derive the per-model verdict
    threshold. Produced by `nfsim-rm`'s `dev/compute_noise_floor.py`.
  - `replicates/{model}/rep_NNN_{model}.gdat` — raw per-rep trajectories.
  - `summary.tsv` — per-model wall-time summary, used by
    `benchmark_full.py` for progress sorting.

## History

| Date | Commit | Scope | Notes |
|---|---|---|---|
| 2026-04-10 | `2de71bf` | All 71 models | First full regen with `-bscb` + locally-reverted upstream #61 bug |
| 2026-04-12 | `05a9a12` | `rm_tlbr_rings` | NFsim reference replaced with hand-rolled SSA (NFsim under-counts on disjoint+symmetric ring patterns) |
| 2026-04-13 | `9cfe36f` | `toy_jim` | NFsim reference replaced with hand-rolled SSA (NFsim over-counts `K(Y~P)` by ~12%) |
