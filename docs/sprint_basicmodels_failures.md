# Sprint outcome — basicmodels parity failures

**Status: complete (2026-04-27).** Started 2026-04-27 against a
26 PASS / 5 FAIL / 2 NO_MATCH baseline on the 33-model basicmodels
parity suite. Closed 2026-04-27 with **29 PASS / 0 FAIL / 0 NO_MATCH
on the trimmed 29-model suite**. Three engine bugs fixed, four
upstream NFsim tests dropped as not-applicable.

## Final state

```
basicmodels @ 5 reps : 29 PASS / 0 FAIL / 0 NO_MATCH (29 imported tests)
ctest                : 1/1
smoke benchmark      : 8/0
feature_coverage     : 55/0  (added 2 regression tests)
corpus guard tier    : 29/0
```

## What was queued at the start

Phase-1 triage (RM @ 5 then @ 20 reps) found these five FAILs:

| model | tz_max @ 5 reps | tz_max @ 20 reps | symptom |
|---|---:|---:|---|
| r02 | 31.96 | 66.06 | `Xpp` over-counted ~1.5×; doubly-phos X accumulates faster than NFsim |
| r07 | 3755.18 | 4039.72 | `full` (5-mol palindrome) under-counted by exactly 2× |
| r18 | 16.74 | 38.97 | `L2(r,r,r)` triply / doubly / singly distribution skew |
| r33 | 177.52 | 177.52 | bond-swap rule plateaus at AC=149; NFsim drives to AC≈5 |
| r35 | 19.47 | 35.97 | `bound_b1` doesn't decay in RM; NFsim drives to 0 |

Plus two NO_MATCH (r31, r34) — no observable to compare against.

## What landed

### Three engine fixes

| commit | subject | model | what it fixed |
|---|---|---|---|
| `5d90724` | `count_multi_molecule_embeddings: branch over all sym partner embeddings` | r07 | Multi-mol Molecules observable BFS committed to the first valid partner embedding and dropped sym-equivalent alternatives. Replaced BFS body with a recursive enumerator that branches over every partner embedding consistent with the walked bond. |
| `7472a07` | `count_multi_mol_fast: forward reacting_local to seed dedup` | r02 | Multi-mol unimolecular rule-rate path called `count_embeddings_single` for the seed without `reacting_local`, so sym embeddings that all map reacting components to the same host slots were never deduped. Threaded `reacting_local` through. |
| `3423b0d` | `init_rule_states: drop redundant compile-time embedding correction` | r18 | The previous commit's seed dedup made `compute_embedding_correction[_multimol]` redundant; with both active the rule rate halved on patterns whose seed had sym non-reacting components. Set `embedding_correction_a/_b = 1.0` for multi-mol (mirrors single-mol) and removed ~110 lines of dead code. |

### Four upstream NFsim tests dropped as not-applicable

| commit | subject | tests |
|---|---|---|
| `85feae1` | `basicmodels: drop r33 and r35 — NFsim-only quirks, BNG2 strict refuses` | r33 (NFsim issue #22 occupied-site bond error); r35 (NFsim issue #14 dot-product split) |
| `9fb2efb` | `basicmodels: drop r31 and r34 — no comparable observables` | r31 (crash test, no observables block); r34 (NFsim issue #24, observable commented out by author) |

For r33 and r35, BNG2.pl `generate_network` was the deciding voice:
on the BNGL the upstream tests carry, BNG2.pl emits the chemistry-correct
network (bound by free-B pool for r33; zero reactions for r35) and
RuleMonkey matches BNG2. The NFsim references they were generated from
captured the historic NFsim quirks they were authored to test, and
including them in the parity suite would have meant verifying RuleMonkey
reproduces NFsim bugs.

The full per-test rationale (with the BNG2-confirmed counter-evidence for
r33 / r35 and the author-comment evidence for r31 / r34) lives in
`tests/reference/basicmodels/PROVENANCE.md`.

### Two regression tests added

Pinned to `tests/models/feature_coverage/`:

- `ft_multimol_sym_obs.bngl` — pre-fix tz=506.91, post-fix tz=3.33.
  Catches the r07 shape (multi-mol Molecules observable on a
  palindromic 5-mol pattern with sym components on the partner).
- `ft_multimol_unimol_unbind_sym.bngl` — pre-fix tz=39.80, post-fix
  tz=2.11. Catches the r02 shape (multi-mol unimolecular unbind rule
  on a host with sym components, where operations don't differentiate
  the embeddings).
- `ft_multimol_pattern_sym_nonreacting.bngl` — pre-fix tz=22.28,
  post-fix tz=1.53. Catches the r18 shape (the embedding-correction
  vs reacting_local-dedup double-correction).

(Three tests landed; a redundant invariant block was tightened during
review on the first one.)

## Suite cleanup

- `tests/reference/basicmodels/PROVENANCE.md` rewritten to lead with
  what the suite *is* (29 imported tests, what they cover, how
  references are generated). The seven upstream NFsim tests not
  carried over (r27, r28, r31, r33, r34, r35, r36) live in an
  appendix grouped by reason — features RM doesn't implement / NFsim
  flags RM doesn't expose / NFsim regressions for bugs that contradict
  BNGL strict / tests with no comparable observables.
- `harness/basicmodels.py` and `harness/generate_basicmodels_refs.py`
  docstrings tightened: a one-line origin pointer to PROVENANCE,
  no exclusion enumeration in the docstring itself.

## Outcomes vs original plan

The original Phase-2 fix order was r07 → r02 → r18, with r33 and r35
queued for after. The cluster hypothesis ("r02 / r07 / r18 might
share a single sym-K root cause") held in spirit — all three lived
in the embedding-counting / seed-dedup machinery — but landed as
three distinct fixes, each on a different code path.

For r33 and r35 the original plan assumed they were RM bugs; both
turned out to be NFsim-specific quirks. Verifying via BNG2.pl
`generate_network` was the bisection step: the moment BNG2 either
generated zero reactions (r35) or generated reactions bound by free-B
count (r33), the divergence with NFsim's reference clearly belonged
to NFsim, not RuleMonkey.

## Re-entry checklist for future basicmodels sprints

1. `cd ~/Code/RuleMonkey && git pull && git status` (clean tree)
2. `cmake --build --preset release && ctest --preset release` (1/1)
3. `python3 harness/benchmark_full.py --tier smoke` (8/0)
4. `python3 harness/basicmodels.py --reps 5` (expect 29/0/0)
5. Open this doc; if a new test is failing, triage first
   (re-run @ 20 reps), then BNGL-bisect, then code.
