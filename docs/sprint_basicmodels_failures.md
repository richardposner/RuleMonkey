# Sprint plan — basicmodels parity failures

**Status: planned, not yet started.** Intended to be opened in a fresh session.

## What landed before this sprint

The basicmodels suite (33 testable models, see `tests/models/nfsim_basicmodels/`)
was rewired to use the same z-score / `tz_max < T_model` verdict as
the corpus harness. References are vendored at
`tests/reference/basicmodels/{xml,ensemble,sim_params.tsv,PROVENANCE.md}`.

First parity run (5 RM reps, no per-model `noise_floor.tsv`, so
T = 5.0 across the board):

```
26 PASS / 5 FAIL / 2 NO_MATCH
```

The 2 NO_MATCH (r31, r34) have no `begin observables` block in their
.bngl source — there's nothing to integrate or compare. Not failures.

The 5 FAIL'd models are this sprint's target.

## The five failures

| Model | tz_max | Worst obs | What the model tests |
|---|---:|---|---|
| r02  | 31.96   | (TBD) | multiple identical phosphorylation sites |
| r07  | 3755.18 | (TBD) | cross-phosphorylation of linked molecules, more complex setup |
| r18  | 16.74   | (TBD) | component handling with basic reactions in molecules with symmetric components |
| r33  | 177.52  | (TBD) | NFSIM ONLY issue22/21 — operations order at low UTL |
| r35  | 19.47   | (TBD) | NFSIM ONLY issue14 — rule disambiguation in bonded complex |

`tz_max` here is a per-rep time-integrated z-score against the
NFsim 100-rep ensemble. T = 5.0 default threshold (from
`benchmark_full.VERDICT_TZ_DEFAULT`); with a calibrated
`noise_floor.tsv` we'd get per-model T = max(5.0, 1.2 × tz_p99).
**r07 at 3755 is two orders of magnitude beyond any plausible
noise floor — that one is structural.**

## Suggested investigation order

I'd start with phase-1 triage (~30 min total) before committing
to a fix order, because some of these may share a root cause and
fixing one could cascade.

### Phase 1: triage

For each of the 5 failures:

1. **Bump RM reps from 5 to 20**, re-run with
   `python3 harness/basicmodels.py r02 r07 r18 r33 r35 --reps 20`.
   Stochastic noise scales as 1/√N; if a model's tz_max drops
   below 5 with more reps, it was a noise artifact at low rep count
   and the model is fine.

2. **Read the .bngl source** (`tests/models/nfsim_basicmodels/v{NN}.bngl`)
   for unusual BNGL features: multi-mol reactant patterns, symmetric
   components, ring-closure rules, unbinding patterns with `+` on the
   product side, large UTL needs, function rate laws.

3. **Spot-check a single observable trajectory**:
   - Run `rm_driver` once with seed=1 on the cached XML
     (`tests/reference/basicmodels/xml/r{NN}.xml`).
   - Compare to `tests/reference/basicmodels/replicates/r{NN}/rep_001.gdat`
     by eye for the worst observable (look for early-time disagreement —
     usually the smoking gun).

4. **Classify** as:
   - `noise` — tz drops below 5 at 20 reps; close the issue
   - `structural` — clear systematic divergence visible in a single rep
   - `flag-related` — divergence depends on a translated NFsim flag
     (revisit `STRIP_FLAGS` in `harness/generate_basicmodels_refs.py`)

### Phase 2: fix in priority order

Suggested order (subject to revision after triage):

1. **r07** first. Its tz=3755 is so far above noise that it's the
   strongest signal of a structural bug. Cross-phosphorylation of
   linked molecules involves multi-mol unimolecular patterns —
   that's a code path with a long bug history (see
   `nfsim-rm` memory entries `project_an_anx_fix.md`,
   `project_multimol_sprint.md`, `project_obs_fm_landed.md`).
   Fixing this bug may incidentally fix r02 / r18 if they share
   the same code path.

2. **r33** next (tz=177.52). The model name calls out "operations
   order and UTL low." Even though we drop `-utl 3` in flag
   translation, the .bngl may have rules whose `max_pattern_size`
   creates an auto-UTL that doesn't match the test's assumptions.
   Worth checking whether RM and NFsim arrive at the same auto-UTL
   for this model.

3. **r02, r18, r35** — if not already fixed by 1 and 2, investigate
   each. r35's "issue14" is rule disambiguation in a bonded complex,
   which intersects with the multi-mol pattern matching in `count_*`
   functions. r18's "symmetric components" intersects with the
   sym-K matching code. Both are in the same neighborhood as r07.

For each fix:

- BNGL bisect first per
  `~/.claude/projects/-Users-wish-Code-RuleMonkey/memory/feedback_bisect_before_reading.md`.
  Strip features until parity is restored — that pinpoints the
  triggering feature without any code reading.
- Reproduce as a focused test in `tests/models/feature_coverage/`
  if a small repro fits there. (One of the corpus's strengths is
  having minimal feature-targeted models; a basicmodels regression
  can become a permanent feature_coverage entry.)
- Land the fix; verify both `harness/benchmark_full.py --tier guard`
  (29/0 baseline) and `harness/basicmodels.py` (this sprint's
  failures should drop) stay green.

### Phase 3: close out

- Update this doc with outcomes per model (passed / fixed / skipped /
  documented as gap).
- Update `docs/FAILING_MODELS.md` with the basicmodels verdicts.
- Bump CHANGELOG with a "Fixed" section listing what landed.
- Consider whether to also generate a calibrated `noise_floor.tsv`
  for the basicmodels suite. Without it the verdict uses T=5
  uniformly; per-model T from self-split calibration would tighten
  signal/noise. ~1 hour of work via porting `nfsim-rm/dev/compute_noise_floor.py`
  to read from `tests/reference/basicmodels/replicates/`.

## Expected outcomes

Realistic ranges:

- **2-3 models PASS at 20 reps** (noise artifacts at 5 reps).
- **1-2 models become real RM bug fixes** (probably in the multi-mol
  matching code path; will land as separate commits with
  feature_coverage regression tests).
- **0-1 models documented as gaps** (e.g., if a basicmodels test
  exercises a feature genuinely outside RM's design — though this
  is unlikely after the r27/r28/r36 cleanup).
- **0 models lost to harness wiring** (the harness is now solid;
  if it breaks, that's a separate bug).

A "good" sprint outcome is **30+/0 PASS** on the next run, with
the residual 0-3 documented in `FAILING_MODELS.md` and the gap
pattern noted in CHANGELOG.

## Re-entry checklist for the fresh session

1. `cd ~/Code/RuleMonkey && git pull && git status` (clean tree)
2. `cmake --build --preset release && ctest --preset release` (1/1)
3. `python3 harness/benchmark_full.py --tier smoke` (8/0 baseline)
4. `python3 harness/basicmodels.py r02 r07 r18 r33 r35 --reps 20` (triage)
5. Open this doc and `docs/FAILING_MODELS.md`; pick a target.

Memory pointer: `project_basicmodels_failures.md` in the new RuleMonkey
auto-memory dir.
