# `edg_*` stress models — rationale and probe targets

These models were authored to break RM (honestly).  Each targets a
specific code path or feature combination that prior coverage left
under-tested, with one or more failure modes the author had reasonable
doubts about.  All run to completion in <2 s under both BNG2 ODE and
RM, with small populations and short t_end.

## Authoring constraints

- Each model includes a `simulate({method=>"nf",...})` action so the
  existing `harness/benchmark_feature_coverage.py` runner picks it up.
- Each carries `# invariant:` declarations where conservation /
  balance laws apply (the runner cross-checks).
- Populations ≤ 60 and t_end ≤ 50 keep total wall under 5 s for the
  full set, even with 20+ NFsim reps.
- BNG2 2.9.3 quirk: `exclude_reactants(N,Pat)` does **not** tolerate a
  space after the comma in this version — keep the form
  `exclude_reactants(1,A(...))`.

## Reference-data policy

Three options, listed by preference for these models:

1. **ODE** (`generate_network` + `simulate_ode`) — deterministic, exact
   in the mean-field limit.  This is the gold standard for the 13
   tightly-bounded models below.
2. **SSA** (`generate_network` + `simulate_ssa`, averaged over ~10
   reps) — matches RM's stochastic semantics including variance.
3. **NFsim** (averaged over ~20 reps) — only useful for cross-check,
   carries z-score noise.

Polymer-ish models (`edg_pattern_local_fcn`, `edg_ring_break_constraint`,
`edg_seeded_ring`) generate larger networks; ODE still works but may
need higher `max_iter` / `max_agg` flags on `generate_network`.

## Model index

| # | Model | Probe target |
|---|---|---|
| 1 | edg_state_increment_chain | Non-alphabetic state declarations (`A(s~s2~s0~s4~s1~s3)`); 8 transitions exercise the canonical-sort path at `simulator.cpp:559-575`. Code that confuses declaration order with sorted index lands on the wrong state. |
| 2 | edg_synth_bonded_complex | `0 -> A(b!1).B(a!1)` — AddBond between *two* newly-added molecules, plus a multi-molecule `DeleteMolecules` for the dimer. Dual `product_mol_a`/`product_mol_b` resolution path has no fast-path coverage today. |
| 3 | edg_time_dependent_rate | Rate uses `time()` (`prod_rate = k0 * exp(-k_decay * time())`). Catches stale-time bugs where the time variable isn't refreshed per propensity recompute. |
| 4 | edg_pattern_local_fcn | Local function with the tag at the **pattern** level (`%x:A(r!1).A(l!1)`). Forces the complex-wide branch of `local_arg_is_molecule=false` (`simulator.cpp:1188-1199`); `ft_local_functions` only exercises the molecule-tag form. |
| 5 | edg_homotrimer_binding | `A(b,b,b)` — three-fold homotypic component. Triggers `has_symmetric_components` and degenerate `count_embeddings_single` past the 2-symmetric case. |
| 6 | edg_ring_break_constraint | Coexisting 2-mer (clean break) and 4-mer chain (still-connected after break) populations under strict `+` semantics — must accept the dimer fire and refuse the 4-mer fire. |
| 7 | edg_compound_op_swap | One rule with **five** distinct OpTypes (StateChange + AddBond + DeleteBond + AddMolecule + DeleteMolecule) in a single `ListOfOperations`. Densest single-rule action profile in the suite. |
| 8 | edg_self_dimerize | `A(s)+A(s) -> A(s!1).A(s!1)` with an **asymmetric** reverse rule that targets one of the two embeddings of the dimer — known trip-hazard for partial-symmetric pattern dedup (msite/fceri history). |
| 9 | edg_state_wildcard_set | `~?` on the **rule** LHS (not just observables) combined with `exclude_reactants`. Distinct code path through the count_embeddings state filter. |
| 10 | edg_multi_pattern_obs | Observables with five alternative patterns each, including a strict-subset overlap (`A(s!+)` and `A(s!+,t!+)`) so doubly-bound A counts twice in Molecules-type but once in Species-type. Tests both summation semantics. |
| 11 | edg_deep_param_chain | Five-deep parameter chain (`rate -> a5 -> a3 -> a4 -> {a1,a2,base}`) declared in non-topological order. BNG2 promotes them to a function chain; RM must recursively resolve. `ft_nested_functions` only goes 3 deep. |
| 12 | edg_seeded_ring | Pre-formed 3-cycle ring (`A(l!1,r!2).A(l!2,r!3).A(l!3,r!1)`) in seed species — every molecule has both bonds; no loose end for the species-init bond resolver. |
| 13 | edg_oscillator | Repressilator with Hill-form rates (n=4) — three coupled negative feedbacks force the rate-recompute path on every event; rate-staleness shows up as dampened oscillations. |
| 14 | edg_fixed_competition | `$E` consumed by three competing rules + a second Fixed type (`$Trash`) — exercises multi-rule replenishment and the "multiple Fixed of distinct types" allowed path. |
| 15 | edg_three_mol_pattern | Single rule with a 3-molecule LHS (`A.B.C`) — ternary embedding. Coverage tops out at 2-mol patterns elsewhere. |
| 16 | edg_dynamic_rate_zero_obs | Rate function with `if(Pool>0, ..., 0)` and a divide-by-Pool form. Catches bugs where a stale "0 propensity" rule never wakes up after the gating observable transitions through 0. |
| 17 | edg_double_state_change | Single rule that flips TWO components on the same molecule simultaneously. Two StateChange ops in one ListOfOperations on distinct components — a rare combination. |
| 18 | edg_branched_polymer | `A(c,c,c,c)` tetravalent component growing branched aggregates with one bond rule. Exercises the cycle-bond / agent-pool BFS shortcut on a topology with degree-4 internal nodes. |
| 19 | edg_zero_rate_rule | A reversible rule whose forward rate is parameter-bound to 0 — propensity tracking must not blow up; the priority queue must skip the rule cleanly. |
| 20 | edg_synth_bind_existing | `A(b) + 0 -> A(b!1).B(a!1)` — synthesise a new B AND bond it to an existing A.  AddBond between an EXISTING reactant and a NEWLY-ADDED product molecule, the third combination on top of the two covered by `edg_synth_bonded_complex` (both new) and ordinary AddBond (both existing). |

## Per-model invariants worth checking

`benchmark_feature_coverage.py` parses `# invariant:` lines.  Models
without explicit invariants rely on the RM-vs-reference verdict only.
The following models have invariants that should hold *every step*:

- conserved totals: 1, 2, 4, 5, 7, 8, 14, 15, 17, 18, 19, 20
- balance: (none in this set yet — candidates for a future revision)

## Verdict reference (verified — `reference/ode/edg_*.gdat` is in tree)

**18 / 20 — ODE gold standard.**  References generated via
`generate_network` + `simulate_ode` and written as TSV to
`reference/ode/edg_*.gdat`.  Deterministic, exact in the mean-field
limit, no stochastic noise.

**2 / 20 — 100-rep NFsim gold standard.**  `edg_pattern_local_fcn` and
`edg_branched_polymer` are network-free by design (chain growth,
tetravalent branching).  Truncated `generate_network` was tried and
rejected — RM doesn't truncate, so a truncated ODE is a different
model.  References live at
`reference/nfsim/edg_pattern_local_fcn.{mean,std}.tsv` and
`reference/nfsim/edg_branched_polymer.{mean,std}.tsv`, each averaged
over 100 NFsim reps with seeds 1..100 (parallel, ~0.5 s wall on this
host).  RM-vs-NFsim verdict uses z-score against the recorded mean ± std.

20-rep RM vs ODE rel-err on the worst observable at t_end (sanity
spot-check; populations small so most are limited by stochastic noise,
not bugs):

| Model | worst-obs rel-err | notes |
|---|---|---|
| edg_state_increment_chain | 0.55 | 5-state tail, low pop; converges with more reps |
| edg_synth_bonded_complex | 0.03 | clean |
| edg_time_dependent_rate | 0.09 | clean (time()-driven decay reference matches) |
| edg_pattern_local_fcn | — | NFsim ref (100 reps). 500-rep investigation found a REAL RM bug: A+A self-binding under bscb under-binds by ~1.5% vs ODE (NFsim ~0.9% bias in same direction). FIXED 2026-04-28 at engine.cpp:2731-2755. After fix, all 4 simulators (BNG2 ODE, BNG2 SSA, NFsim, RM) agree within stochastic noise on every edg_* model. See `project_self_binding_bias.md` in auto-memory. Same bug also affected existing corpus model `combo_addbond_connected`. |
| edg_homotrimer_binding | 0.64 | A_act, 20-molecule stochastic — re-check at 100 reps |
| edg_ring_break_constraint | 0.07 | clean |
| edg_compound_op_swap | 0.00 | clean |
| edg_self_dimerize | 0.46 | A_free, ~8/40 at equilibrium — Poisson noise dominates |
| edg_state_wildcard_set | 0.02 | clean |
| edg_multi_pattern_obs | 0.03 | clean |
| edg_deep_param_chain | 0.00 | clean — RM resolves the function chain correctly |
| edg_seeded_ring | 0.00 | clean (with the spin rule added) |
| edg_oscillator | 0.09 | clean for an oscillating system at small N |
| edg_fixed_competition | 0.00 | clean — Fixed clamp holds |
| edg_three_mol_pattern | 0.06 | clean |
| edg_dynamic_rate_zero_obs | 0.42 | P_n ≈ 0.5 — high relative noise on tiny mean |
| edg_double_state_change | 0.01 | clean |
| edg_branched_polymer | — | NFsim ref (100 reps); RM matches cleanly, z<0.5 on all observables |
| edg_zero_rate_rule | 0.00 | clean — zero-rate reverse correctly inert |
| edg_synth_bind_existing | 0.18 | within noise |

## Gold-standard generation notes

- `generate_network({overwrite=>1})` worked out of the box for 17/20.
- `edg_seeded_ring` initially produced 0 reactions (rule didn't match
  any seed), so a spin-flip rule was added for non-trivial dynamics
  while preserving the seed-species topology test.
- The 2 polymer-style models (`edg_pattern_local_fcn`,
  `edg_branched_polymer`) are network-free in nature.  Truncated
  `generate_network` produces a different model than RM simulates —
  rejected as the gold standard.  These two use 100-rep NFsim refs
  instead (`reference/nfsim/<model>.{mean,std}.tsv`).
- BNG2 2.9.3 quirk found: `exclude_reactants(N,Pat)` does not tolerate
  a space after the comma.  Use `exclude_reactants(1,A(...))`.
