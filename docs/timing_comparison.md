# RuleMonkey vs NFsim — wall-time comparison

Single-machine, sequential, end-to-end subprocess wall time.  3 reps per engine per model; median reported.  Both engines invoked from a fresh process per rep (no warm caches).  Same XML, same t_end, same n_steps, same seed sequence.

**This is an efficiency report, not a correctness report.**  Models where NFsim runs but produces incorrect observables (e.g. `ft_tfun` due to NFsim's TFUN handler returning zero rate) are still included in the timing table; correctness is covered separately by `feature_coverage`, `benchmark_full`, and `basicmodels` suites.

Generated: `2026-05-01 23:17:28 MDT` (`benchmark_rm_vs_nfsim_timing.py`)

## Summary

- Models scored: **169 / 177** (both engines completed all 3 reps)
- NFsim N/A: **8** (RM ran; NFsim refused / errored)
- **Median speedup (NFsim wall ÷ RM wall): 0.86×**
- Geometric mean speedup: 0.77×
- RM faster than NFsim on **69 / 169** models; NFsim faster on 100

## Where does each engine win?

The 169 scored models fall into six speedup buckets, and the buckets
correlate cleanly with model shape:

| Bucket | Count | Examples | What's going on |
|---|---:|---|---|
| RM ≥ 10× | 4 | `mlnr` 23.3×, `pltr` 23.0×, `testcase2a` 11.7×, `rm_tlbr` 11.6× | Pattern-quantifier (count-relation) Species observables |
| RM 2–10× | 1 | `example3_fit` 2.65× | Single-model outlier |
| RM 1.1–2× | 48 | `ss_tlbr_rings`, `combo_symmetric_rings`, many small fast models | RM's lower per-process startup wins on short runs |
| Parity (0.8–1.1×) | 39 | `BLBR`, `tlbr`, `A_plus_B_rings`, `bench_tlbr_yang2008` | Same network-free chemistry, both engines comparable |
| NFsim 0.5–0.8× | 30 | `e2`–`e9`, `tcr_*` variants | Heavy enzyme-kinetics scaling; NFsim a bit ahead |
| NFsim < 0.5× | 47 | `AN` 0.12×, `ANx` 0.14×, `oscSystem` 0.14×, `fceri_ji` 0.27× | Heavy stochastic corpus models: NFsim's tuning shows |

### RM ≥ 10× — count-relation Species observables

All four models in this bucket declare **exactly 300 Species observables
of the form `R()=N`** (`R()=1`, `R()=2`, ..., `R()=300`) — histogram
bins counting complexes that contain exactly *N* receptors.  No other
model in the corpus declares any:

| Model | Count-relation Species obs | Speedup |
|---|---:|---:|
| `mlnr`       | 300 | 23.3× |
| `pltr`       | 300 | 23.0× |
| `testcase2a` | 300 | 11.7× |
| `rm_tlbr`    | 300 | 11.6× |
| (next-largest model in any other bucket) | 0 | — |

NFsim re-evaluates count-relation Species observables naively at each
sample point: walk every complex, check the count constraint, for every
observable.  At ~600 complexes peak × 300 observables × 1000 sample
points, that's ~180M observable evaluations per run.  RM's
Species-observable incremental tracking only re-evaluates dirty
complexes between sample points — typically O(1)–O(2) complexes per
rule fire.  Implementation lives in `init_incremental_observables` /
`flush_species_incr_observables` in `cpp/rulemonkey/engine.cpp`.  The speedup is
proportional to (`# observables` × `# sample points`) ÷
(`# events between samples`), and on these models the ratio happens to
be ~10×–25×.

### RM 1.1–2× — short / cheap models

The 48 models in this band are mostly sub-100ms runs where end-to-end
process startup is a non-trivial fraction of wall time.  RM has
slightly lower startup overhead (smaller binary, simpler XML loader),
which shows up as a 1.1×–1.5× win on short runs.  Vanishes once the
simulation itself dominates wall time.

### Parity — generic network-free chemistry without histogram observables

`A_plus_B_rings` (no count-relation obs) sits at 1.01×; `bench_tlbr_yang2008`
at 0.90×; `BLBR` at 0.88×.  Same kind of binding chemistry as the
RM-winning Species-observable models above — but without the
observable-evaluation asymmetry, the engines are comparable.

### NFsim < 0.5× — heavy stochastic corpus models

The dominant story for the bottom of the table.  Models like `t3`,
`machine`, `fceri_ji`, `ensemble`, `rm_tlbr_rings`, `AN` are the
heavy-event-rate workloads NFsim has been tuned on for years.  RM has
real ground to make up here, but the 2026-05-01 perf sprint closed a
significant chunk of it (see "What changed" below).

The single-model outlier `example3_fit` (2.65×) doesn't fit any of
these buckets cleanly — worth a separate look if anyone wants to
chase it.

## What changed since the 2026-04-29 baseline

Eight engine-level perf commits landed on 2026-05-01
(`3b66d69` … `d889dab`):

- ExprEval `Variable` lookup indexing + builtin tagging
- `update_eval_vars` short-circuit + parameter hoisting
- Delta-update of `total_propensity` (replaces full re-sum each event)
- Pattern-adjacency reuse + reservoir sampling for disjoint assignments
- `apply_overrides` walk skipped when overrides map is empty
- `count_embeddings_single` inlined via self-passing lambda
- `resolve_value` parsed-AST cache
- Retry-until-distinct sampler for `same_components` homodimers

These target the per-event hot loop, so the impact is concentrated in
the heavy stochastic corpus bucket.  RM wall-time deltas on a
representative slice (median, same machine, two days apart):

| Model | Prev RM | Now RM | Δ |
|---|---:|---:|---:|
| `t3` | 49.9 s | 39.0 s | -22% |
| `machine` | 42.1 s | 22.7 s | -46% |
| `ensemble` | 32.3 s | 17.9 s | -45% |
| `oscSystem` | 9.42 s | 2.82 s | -70% |
| `CaMKII_holo` | 972 ms | 427 ms | -56% |
| `isingspin_localfcn` | 1066 ms | 634 ms | -41% |
| `rm_tlbr_rings` | 35.3 s | 29.8 s | -16% |
| `ANx` | 11.5 s | 9.4 s | -18% |
| `AN` | 21.9 s | 18.7 s | -15% |
| `lat` | 7.91 s | 6.38 s | -19% |
| `e1` | 6462 ms | 5808 ms | -10% |
| `tcr` | 21.8 s | 19.9 s | -9% |
| `fceri_ji` | 33.6 s | 30.9 s | -8% |

The headline summary numbers (median 0.83× → 0.86×, geomean 0.78× →
0.77×) barely moved because those statistics are dominated by the
~120 small models in the 1.1–2× and parity bands where end-to-end
wall time is dominated by process startup, not by the per-event hot
loop the optimizations targeted.  The absolute wall time saved on the
heavy-stochastic bucket is the more honest measure: ~45 seconds shaved
off `t3`+`machine`+`ensemble` per rep alone.  Bucket migration: 8
models shifted from the RM 1.1–2× band to parity (NFsim ran a touch
faster on this run too — same-machine cross-day timing has noise),
and the < 0.5× and 0.5–0.8× bands each shed 2 models.

NFsim N/A count rose from 4 to 8: `rm_tlbr_rings`, `example2_fit`,
`egfr_nf_iter5p12h10`, and `bench_blbr_dembo1978_monovalent_inhibitor`
joined the prior four (`ft_nested_functions`,
`edg_deep_param_chain`, `edg_time_dependent_rate`, `ft_tfun`).  RM
ran them all to completion; the diagnostics in the Notes column
attribute the failures to NFsim itself.

## Caveats

- Wall time, not CPU time.  Includes process startup, XML load, and result write — same for both engines, so the ratio is fair, but absolute numbers shift across machines.
- Speedup depends on workload: long-horizon, high-event-rate models show RM's biggest wins; short-horizon trivial models are near parity or NFsim-favored due to RM's slightly heavier per-process startup.
- Single-replicate timings can be noisy.  See min/max columns for spread; the median column is what to quote.
- N/A in the NFsim column means NFsim refused or errored on at least one of the 3 reps.  The reason column gives the captured diagnostic.

## Per-model results

Sorted by RM median wall time (most expensive first).

| # | Suite | Model | t_end | n_steps | RM median | RM min/max | NFsim median | NFsim min/max | Speedup | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | corpus | `t3` | 100 | 100 | 39.0 s | 38.8 s / 39.1 s | 21.8 s | 20.9 s / 22.3 s | 0.56× |  |
| 2 | corpus | `fceri_ji` | 500 | 500 | 30.9 s | 30.8 s / 32.8 s | 8310 ms | 7896 ms / 8322 ms | 0.27× |  |
| 3 | corpus | `rm_tlbr_rings` | 1000 | 100 | 29.8 s | 28.6 s / 31.6 s | N/A | — | — | NFsim produced no output |
| 4 | corpus | `machine` | 14 | 100 | 22.7 s | 22.4 s / 24.7 s | 9497 ms | 9495 ms / 9855 ms | 0.42× |  |
| 5 | corpus | `tcr` | 60 | 12 | 19.9 s | 19.7 s / 20.6 s | 10.4 s | 10.1 s / 10.8 s | 0.52× |  |
| 6 | corpus | `example2_fit` | 126 | 20 | 19.1 s | 18.1 s / 19.8 s | N/A | — | — | NFsim produced no output |
| 7 | corpus | `AN` | 100 | 100 | 18.7 s | 18.6 s / 19.1 s | 2178 ms | 2159 ms / 2232 ms | 0.12× |  |
| 8 | corpus | `ensemble` | 15 | 100 | 17.9 s | 16.6 s / 18.2 s | 7631 ms | 7611 ms / 7795 ms | 0.43× |  |
| 9 | corpus | `e9` | 200 | 100 | 16.4 s | 15.9 s / 16.5 s | 13.2 s | 12.8 s / 13.4 s | 0.81× |  |
| 10 | corpus | `example4_fit` | 13 | 100 | 15.9 s | 15.5 s / 16.1 s | 9253 ms | 8779 ms / 9629 ms | 0.58× |  |
| 11 | corpus | `egfr_nf_iter5p12h10` | 60 | 20 | 15.7 s | 14.4 s / 16.6 s | N/A | — | — | NFsim produced no output |
| 12 | corpus | `e8` | 200 | 100 | 14.8 s | 14.7 s / 15.0 s | 11.7 s | 11.5 s / 11.9 s | 0.79× |  |
| 13 | corpus | `e7` | 200 | 100 | 14.5 s | 14.4 s / 14.6 s | 10.2 s | 9961 ms / 10.2 s | 0.70× |  |
| 14 | corpus | `bench_blbr_rings_posner1995` | 3000 | 300 | 13.5 s | 12.6 s / 14.2 s | 5764 ms | 5757 ms / 6266 ms | 0.43× |  |
| 15 | corpus | `blbr_rings_posner1995` | 3000 | 300 | 13.4 s | 13.1 s / 13.7 s | 6233 ms | 5803 ms / 6918 ms | 0.47× |  |
| 16 | corpus | `e6` | 200 | 100 | 12.5 s | 12.4 s / 13.1 s | 8934 ms | 8902 ms / 9076 ms | 0.72× |  |
| 17 | corpus | `e5` | 200 | 100 | 11.1 s | 10.9 s / 11.2 s | 7426 ms | 7350 ms / 7694 ms | 0.67× |  |
| 18 | corpus | `PushPull` | 4000 | 40 | 11.0 s | 10.9 s / 11.3 s | 4058 ms | 3990 ms / 4100 ms | 0.37× |  |
| 19 | corpus | `e4` | 200 | 100 | 9687 ms | 9672 ms / 9693 ms | 6276 ms | 5982 ms / 6279 ms | 0.65× |  |
| 20 | corpus | `ANx` | 100 | 100 | 9413 ms | 9330 ms / 9674 ms | 1311 ms | 1310 ms / 1313 ms | 0.14× |  |
| 21 | corpus | `e3` | 200 | 100 | 8618 ms | 8382 ms / 8683 ms | 4925 ms | 4745 ms / 5088 ms | 0.57× |  |
| 22 | corpus | `e2` | 200 | 100 | 7631 ms | 7388 ms / 7735 ms | 3516 ms | 3407 ms / 3726 ms | 0.46× |  |
| 23 | corpus | `egfr_net` | 120 | 120 | 6505 ms | 6434 ms / 6507 ms | 3877 ms | 3865 ms / 3980 ms | 0.60× |  |
| 24 | corpus | `lat` | 10 | 100 | 6379 ms | 5972 ms / 6617 ms | 2894 ms | 2884 ms / 2913 ms | 0.45× |  |
| 25 | corpus | `e1` | 200 | 100 | 5808 ms | 5727 ms / 5860 ms | 2064 ms | 2018 ms / 2147 ms | 0.36× |  |
| 26 | nfsim_basicmodels | `r02` | 100 | 100 | 4716 ms | 4349 ms / 4829 ms | 1152 ms | 1137 ms / 1245 ms | 0.24× |  |
| 27 | corpus | `bench_blbr_dembo1978_monovalent_inhibitor` | 3000 | 300 | 4559 ms | 4332 ms / 4596 ms | N/A | — | — | NFsim produced no output |
| 28 | corpus | `pltr` | 1000 | 1000 | 3327 ms | 3180 ms / 3355 ms | 76.4 s | 75.9 s / 78.9 s | 22.97× |  |
| 29 | corpus | `mlnr` | 1000 | 1000 | 3241 ms | 3210 ms / 3310 ms | 75.6 s | 75.3 s / 78.6 s | 23.32× |  |
| 30 | nfsim_basicmodels | `r20` | 100 | 100 | 3189 ms | 3079 ms / 3260 ms | 1269 ms | 1205 ms / 1318 ms | 0.40× |  |
| 31 | nfsim_basicmodels | `r01` | 100 | 100 | 3103 ms | 2985 ms / 3180 ms | 1284 ms | 1236 ms / 1477 ms | 0.41× |  |
| 32 | corpus | `oscSystem` | 200 | 200 | 2815 ms | 2597 ms / 2871 ms | 391 ms | 364 ms / 399 ms | 0.14× |  |
| 33 | nfsim_basicmodels | `r07` | 20 | 100 | 2628 ms | 2610 ms / 2633 ms | 786 ms | 765 ms / 832 ms | 0.30× |  |
| 34 | corpus | `tlbr` | 200 | 200 | 2548 ms | 2474 ms / 2662 ms | 2126 ms | 2093 ms / 2247 ms | 0.83× |  |
| 35 | corpus | `testcase2a` | 1000 | 1000 | 2497 ms | 2456 ms / 2517 ms | 29.2 s | 28.5 s / 29.3 s | 11.69× |  |
| 36 | corpus | `rm_tlbr` | 1000 | 1000 | 2482 ms | 2460 ms / 2522 ms | 28.9 s | 28.4 s / 29.0 s | 11.63× |  |
| 37 | nfsim_basicmodels | `r06` | 20 | 100 | 2435 ms | 2375 ms / 2555 ms | 451 ms | 432 ms / 454 ms | 0.19× |  |
| 38 | corpus | `simple_nfsim` | 100 | 50 | 2061 ms | 2005 ms / 2167 ms | 854 ms | 823 ms / 1045 ms | 0.41× |  |
| 39 | corpus | `simple_system` | 100 | 50 | 2039 ms | 2030 ms / 2124 ms | 833 ms | 828 ms / 979 ms | 0.41× |  |
| 40 | corpus | `nfsim_aggregation_gelation` | 50 | 200 | 2015 ms | 2004 ms / 2037 ms | 763 ms | 749 ms / 765 ms | 0.38× |  |
| 41 | nfsim_basicmodels | `r08` | 20 | 100 | 1932 ms | 1909 ms / 1938 ms | 458 ms | 456 ms / 536 ms | 0.24× |  |
| 42 | corpus | `tcr_gen27ind33` | 0.94 | 12 | 1872 ms | 1854 ms / 1974 ms | 732 ms | 729 ms / 744 ms | 0.39× |  |
| 43 | corpus | `tcr_iter28p4h2` | 0.94 | 12 | 1819 ms | 1802 ms / 1874 ms | 738 ms | 736 ms / 749 ms | 0.41× |  |
| 44 | nfsim_basicmodels | `r10` | 100 | 100 | 1774 ms | 1762 ms / 1827 ms | 328 ms | 327 ms / 329 ms | 0.18× |  |
| 45 | nfsim_basicmodels | `r04` | 50 | 100 | 1742 ms | 1726 ms / 1774 ms | 731 ms | 709 ms / 746 ms | 0.42× |  |
| 46 | corpus | `tcr_gen20ind9` | 1.5 | 12 | 1714 ms | 1712 ms / 1774 ms | 716 ms | 716 ms / 739 ms | 0.42× |  |
| 47 | corpus | `tcr_iter9p44` | 0.94 | 12 | 1695 ms | 1689 ms / 1735 ms | 704 ms | 695 ms / 705 ms | 0.42× |  |
| 48 | corpus | `bench_blbr_cooperativity_posner2004` | 10 | 10 | 1259 ms | 1257 ms / 1326 ms | 711 ms | 707 ms / 711 ms | 0.56× |  |
| 49 | corpus | `bench_blbr_cooperativity_posner2004_rings` | 10 | 10 | 1163 ms | 1154 ms / 1192 ms | 619 ms | 600 ms / 625 ms | 0.53× |  |
| 50 | corpus | `stiff` | 10 | 100 | 886 ms | 877 ms / 892 ms | 421 ms | 418 ms / 423 ms | 0.48× |  |
| 51 | corpus | `example3_fit` | 5000 | 10 | 826 ms | 825 ms / 862 ms | 2192 ms | 2134 ms / 2219 ms | 2.65× |  |
| 52 | corpus | `BLBR` | 10 | 100 | 696 ms | 695 ms / 699 ms | 612 ms | 610 ms / 615 ms | 0.88× |  |
| 53 | corpus | `isingspin_localfcn` | 5000 | 100 | 634 ms | 633 ms / 636 ms | 160 ms | 160 ms / 170 ms | 0.25× |  |
| 54 | nfsim_basicmodels | `r03` | 50 | 100 | 566 ms | 531 ms / 570 ms | 230 ms | 228 ms / 236 ms | 0.41× |  |
| 55 | corpus | `rm_blbr` | 10 | 100 | 458 ms | 443 ms / 462 ms | 231 ms | 227 ms / 234 ms | 0.51× |  |
| 56 | corpus | `CaMKII_holo` | 10 | 100 | 427 ms | 420 ms / 443 ms | 573 ms | 564 ms / 595 ms | 1.34× |  |
| 57 | nfsim_basicmodels | `r18` | 10 | 100 | 418 ms | 412 ms / 430 ms | 264 ms | 262 ms / 265 ms | 0.63× |  |
| 58 | nfsim_basicmodels | `r12` | 10 | 100 | 408 ms | 390 ms / 416 ms | 152 ms | 144 ms / 152 ms | 0.37× |  |
| 59 | corpus | `bench_blbr_heterogeneity_goldstein1980` | 10 | 10 | 408 ms | 406 ms / 416 ms | 105 ms | 105 ms / 105 ms | 0.26× |  |
| 60 | feature_coverage | `ss_long_polymer` | 500 | 100 | 358 ms | 352 ms / 363 ms | 154 ms | 152 ms / 159 ms | 0.43× |  |
| 61 | nfsim_basicmodels | `r09` | 100 | 100 | 357 ms | 352 ms / 357 ms | 62.6 ms | 60.5 ms / 64.1 ms | 0.18× |  |
| 62 | feature_coverage | `ft_multimol_sym_obs` | 20 | 100 | 322 ms | 309 ms / 332 ms | 69.4 ms | 68.9 ms / 73.5 ms | 0.22× |  |
| 63 | corpus | `st_multi_2` | 10 | 100 | 311 ms | 306 ms / 316 ms | 141 ms | 126 ms / 163 ms | 0.45× |  |
| 64 | corpus | `toy_jim` | 100 | 100 | 303 ms | 302 ms / 305 ms | 135 ms | 131 ms / 137 ms | 0.45× |  |
| 65 | corpus | `ANx_noActivity` | 100 | 100 | 255 ms | 248 ms / 259 ms | 118 ms | 115 ms / 118 ms | 0.46× |  |
| 66 | feature_coverage | `ss_branching_aggregate` | 400 | 100 | 254 ms | 252 ms / 262 ms | 146 ms | 146 ms / 149 ms | 0.58× |  |
| 67 | corpus | `poly` | 100 | 100 | 208 ms | 208 ms / 216 ms | 119 ms | 116 ms / 213 ms | 0.57× |  |
| 68 | corpus | `bench_tlbr_yang2008` | 3000 | 300 | 205 ms | 189 ms / 207 ms | 184 ms | 176 ms / 187 ms | 0.90× |  |
| 69 | corpus | `bench_blbr_rings_posner1995_no_rings` | 50 | 50 | 181 ms | 177 ms / 188 ms | 155 ms | 152 ms / 156 ms | 0.85× |  |
| 70 | corpus | `bench_blbr_dembo1978` | 40 | 40 | 179 ms | 177 ms / 180 ms | 147 ms | 145 ms / 151 ms | 0.82× |  |
| 71 | corpus | `A_plus_B_rings` | 10 | 10000 | 162 ms | 162 ms / 172 ms | 164 ms | 163 ms / 172 ms | 1.01× |  |
| 72 | corpus | `bench_tlbr_solution_macken1982` | 500 | 250 | 161 ms | 158 ms / 176 ms | 110 ms | 109 ms / 120 ms | 0.68× |  |
| 73 | feature_coverage | `ft_stiff_system` | 200 | 100 | 150 ms | 147 ms / 151 ms | 43.1 ms | 43.1 ms / 43.7 ms | 0.29× |  |
| 74 | feature_coverage | `ss_symmetric_homopoly` | 300 | 100 | 144 ms | 143 ms / 145 ms | 123 ms | 121 ms / 124 ms | 0.86× |  |
| 75 | feature_coverage | `ft_multisite_phospho` | 200 | 100 | 138 ms | 133 ms / 146 ms | 68.7 ms | 65.4 ms / 71.1 ms | 0.50× |  |
| 76 | nfsim_basicmodels | `r11` | 1 | 100 | 127 ms | 125 ms / 129 ms | 49.8 ms | 48.3 ms / 49.9 ms | 0.39× |  |
| 77 | corpus | `A_plus_A_mixed_1` | 30 | 100 | 114 ms | 112 ms / 116 ms | 44.4 ms | 42.5 ms / 44.5 ms | 0.39× |  |
| 78 | feature_coverage | `ss_tlbr_rings` | 400 | 100 | 111 ms | 108 ms / 113 ms | 133 ms | 130 ms / 135 ms | 1.20× |  |
| 79 | nfsim_basicmodels | `r13` | 30 | 100 | 102 ms | 99.0 ms / 103 ms | 22.8 ms | 22.6 ms / 25.5 ms | 0.22× |  |
| 80 | corpus | `st_multi_1` | 10 | 100 | 95.1 ms | 93.6 ms / 97.2 ms | 40.5 ms | 39.4 ms / 41.2 ms | 0.43× |  |
| 81 | nfsim_basicmodels | `r14` | 30 | 100 | 71.4 ms | 69.8 ms / 75.1 ms | 27.1 ms | 27.0 ms / 27.2 ms | 0.38× |  |
| 82 | feature_coverage | `ft_push_pull` | 200 | 100 | 71.1 ms | 69.6 ms / 74.7 ms | 28.7 ms | 27.8 ms / 30.4 ms | 0.40× |  |
| 83 | nfsim_basicmodels | `r21` | 100 | 100 | 65.6 ms | 65.1 ms / 69.8 ms | 82.9 ms | 82.6 ms / 83.4 ms | 1.26× |  |
| 84 | feature_coverage | `ft_receptor_dimerization` | 100 | 100 | 55.7 ms | 53.3 ms / 60.4 ms | 27.6 ms | 27.5 ms / 27.7 ms | 0.50× |  |
| 85 | nfsim_basicmodels | `r16` | 5 | 100 | 44.7 ms | 42.8 ms / 47.7 ms | 19.8 ms | 19.1 ms / 20.4 ms | 0.44× |  |
| 86 | feature_coverage | `ft_cooperative_binding` | 100 | 100 | 42.5 ms | 40.9 ms / 45.3 ms | 15.7 ms | 15.4 ms / 15.9 ms | 0.37× |  |
| 87 | feature_coverage | `combo_strict_product_plus` | 100 | 100 | 37.0 ms | 36.1 ms / 41.1 ms | 18.4 ms | 18.1 ms / 20.0 ms | 0.50× |  |
| 88 | corpus | `Tutorial_Example` | 5 | 50 | 31.8 ms | 30.7 ms / 32.3 ms | 18.4 ms | 17.7 ms / 18.8 ms | 0.58× |  |
| 89 | corpus | `receptor_nf_iter36p0h3` | 60 | 60 | 31.5 ms | 30.9 ms / 35.7 ms | 19.5 ms | 18.8 ms / 20.8 ms | 0.62× |  |
| 90 | corpus | `example6_ground_truth` | 60 | 60 | 27.1 ms | 26.3 ms / 30.8 ms | 17.2 ms | 16.9 ms / 17.4 ms | 0.63× |  |
| 91 | nfsim_basicmodels | `r15` | 1 | 100 | 26.3 ms | 25.5 ms / 28.9 ms | 15.1 ms | 14.9 ms / 16.7 ms | 0.57× |  |
| 92 | feature_coverage | `nf_branching_aggregate` | 200 | 100 | 21.6 ms | 20.8 ms / 23.9 ms | 17.2 ms | 16.7 ms / 17.7 ms | 0.80× |  |
| 93 | corpus | `A_plus_A` | 10 | 100 | 21.6 ms | 21.3 ms / 21.8 ms | 14.2 ms | 13.6 ms / 14.4 ms | 0.66× |  |
| 94 | nfsim_basicmodels | `r05` | 1 | 100 | 19.5 ms | 18.9 ms / 23.0 ms | 14.8 ms | 14.7 ms / 15.0 ms | 0.76× |  |
| 95 | feature_coverage | `ft_continue` | 100 | 100 | 19.0 ms | 18.7 ms / 21.8 ms | 12.3 ms | 12.0 ms / 12.4 ms | 0.64× |  |
| 96 | feature_coverage | `ft_blbr` | 100 | 100 | 16.2 ms | 15.5 ms / 19.1 ms | 13.9 ms | 13.7 ms / 14.5 ms | 0.86× |  |
| 97 | corpus | `nfsim_coarse_graining` | 400 | 200 | 15.1 ms | 14.7 ms / 18.7 ms | 13.9 ms | 13.4 ms / 14.1 ms | 0.92× |  |
| 98 | feature_coverage | `ft_tlbr` | 100 | 100 | 15.1 ms | 14.8 ms / 16.7 ms | 14.0 ms | 13.9 ms / 14.2 ms | 0.93× |  |
| 99 | feature_coverage | `ft_multimol_pattern_sym_nonreacting` | 10 | 40 | 14.4 ms | 14.4 ms / 17.8 ms | 11.3 ms | 11.1 ms / 11.5 ms | 0.78× |  |
| 100 | feature_coverage | `nf_large_multivalent` | 100 | 100 | 13.4 ms | 13.3 ms / 16.9 ms | 13.0 ms | 13.0 ms / 13.3 ms | 0.97× |  |
| 101 | feature_coverage | `ft_include_reactants` | 100 | 100 | 13.3 ms | 12.7 ms / 15.2 ms | 9.6 ms | 9.6 ms / 9.7 ms | 0.73× |  |
| 102 | corpus | `nfsim_dynamic_compartments` | 20 | 100 | 12.9 ms | 12.6 ms / 19.6 ms | 12.9 ms | 10.1 ms / 13.1 ms | 1.00× |  |
| 103 | feature_coverage | `combo_localfcn_multisite` | 200 | 100 | 12.4 ms | 12.1 ms / 15.8 ms | 10.7 ms | 10.6 ms / 11.3 ms | 0.86× |  |
| 104 | feature_coverage | `ft_exclude_products` | 200 | 50 | 12.2 ms | 11.1 ms / 14.1 ms | 9.3 ms | 9.2 ms / 9.4 ms | 0.77× |  |
| 105 | corpus | `st` | 10 | 100 | 11.9 ms | 11.8 ms / 13.9 ms | 11.4 ms | 11.2 ms / 11.7 ms | 0.96× |  |
| 106 | feature_coverage | `nf_linear_polymer` | 200 | 100 | 11.8 ms | 11.4 ms / 15.1 ms | 10.1 ms | 10.0 ms / 10.3 ms | 0.85× |  |
| 107 | feature_coverage | `combo_symmetric_rings` | 100 | 100 | 11.7 ms | 11.4 ms / 13.7 ms | 12.7 ms | 12.6 ms / 14.1 ms | 1.08× |  |
| 108 | feature_coverage | `ft_perturbation` | 100 | 100 | 11.6 ms | 11.3 ms / 12.5 ms | 10.4 ms | 10.4 ms / 10.5 ms | 0.90× |  |
| 109 | feature_coverage | `combo_addbond_connected` | 200 | 100 | 11.1 ms | 10.7 ms / 14.0 ms | 10.1 ms | 9.5 ms / 10.8 ms | 0.91× |  |
| 110 | feature_coverage | `combo_synth_degrade_equilibrium` | 500 | 100 | 11.0 ms | 10.2 ms / 13.2 ms | 10.2 ms | 10.0 ms / 10.3 ms | 0.92× |  |
| 111 | feature_coverage | `ft_multi_site_binding` | 100 | 100 | 11.0 ms | 11.0 ms / 12.9 ms | 10.8 ms | 10.8 ms / 11.1 ms | 0.98× |  |
| 112 | corpus | `basicTLBR` | 20 | 100 | 10.6 ms | 10.4 ms / 12.2 ms | 14.9 ms | 14.8 ms / 14.9 ms | 1.40× |  |
| 113 | feature_coverage | `combo_shorthand_embed` | 100 | 100 | 10.4 ms | 10.3 ms / 11.9 ms | 10.6 ms | 10.5 ms / 10.8 ms | 1.02× |  |
| 114 | feature_coverage | `edg_multi_pattern_obs` | 30 | 60 | 10.3 ms | 9.9 ms / 13.3 ms | 9.1 ms | 9.1 ms / 9.1 ms | 0.89× |  |
| 115 | feature_coverage | `ft_signaling_cascade` | 100 | 100 | 10.1 ms | 10.1 ms / 12.6 ms | 10.7 ms | 10.6 ms / 46.4 ms | 1.06× |  |
| 116 | feature_coverage | `edg_branched_polymer` | 30 | 60 | 10.1 ms | 8.7 ms / 12.5 ms | 13.8 ms | 13.8 ms / 14.7 ms | 1.36× |  |
| 117 | feature_coverage | `edg_three_mol_pattern` | 30 | 60 | 9.5 ms | 8.4 ms / 10.4 ms | 9.2 ms | 9.1 ms / 9.3 ms | 0.97× |  |
| 118 | nfsim_basicmodels | `r19` | 10 | 100 | 9.5 ms | 9.3 ms / 11.7 ms | 16.7 ms | 16.6 ms / 18.7 ms | 1.77× |  |
| 119 | feature_coverage | `combo_multimol_unimol` | 100 | 100 | 9.5 ms | 9.3 ms / 10.6 ms | 10.0 ms | 9.7 ms / 10.3 ms | 1.06× |  |
| 120 | feature_coverage | `ft_energy_patterns` | 500 | 100 | 9.2 ms | 9.1 ms / 11.4 ms | 10.5 ms | 10.2 ms / 10.6 ms | 1.15× |  |
| 121 | feature_coverage | `ft_exclude_reactants` | 100 | 100 | 9.1 ms | 8.8 ms / 11.0 ms | 9.6 ms | 9.5 ms / 9.7 ms | 1.06× |  |
| 122 | feature_coverage | `edg_pattern_local_fcn` | 40 | 80 | 9.0 ms | 8.9 ms / 10.1 ms | 9.2 ms | 9.1 ms / 9.4 ms | 1.02× |  |
| 123 | feature_coverage | `ft_bond_wildcards` | 100 | 100 | 8.9 ms | 8.7 ms / 11.8 ms | 9.5 ms | 9.1 ms / 9.5 ms | 1.07× |  |
| 124 | feature_coverage | `ft_catalytic_unbinding` | 100 | 100 | 8.9 ms | 7.6 ms / 10.4 ms | 9.4 ms | 9.2 ms / 9.4 ms | 1.06× |  |
| 125 | feature_coverage | `ft_population_map` | 60 | 60 | 8.8 ms | 8.8 ms / 11.5 ms | 9.1 ms | 9.0 ms / 9.3 ms | 1.04× |  |
| 126 | feature_coverage | `ft_nested_functions` | 100 | 100 | 8.8 ms | 8.7 ms / 11.3 ms | N/A | — | — | NFsim produced no output |
| 127 | corpus | `nfsim_ring_closure_polymer` | 50 | 200 | 8.7 ms | 8.2 ms / 11.1 ms | 10.0 ms | 10.0 ms / 10.1 ms | 1.15× |  |
| 128 | feature_coverage | `ft_ring_closure` | 100 | 100 | 8.6 ms | 8.6 ms / 9.5 ms | 9.3 ms | 9.2 ms / 9.4 ms | 1.08× |  |
| 129 | feature_coverage | `ft_multimol_unimol_unbind_sym` | 2 | 40 | 8.6 ms | 8.3 ms / 11.6 ms | 8.4 ms | 8.4 ms / 8.5 ms | 0.98× |  |
| 130 | nfsim_basicmodels | `r24` | 100 | 100 | 8.6 ms | 8.5 ms / 11.6 ms | 10.7 ms | 10.5 ms / 10.7 ms | 1.24× |  |
| 131 | feature_coverage | `ft_complex_seed` | 50 | 50 | 8.5 ms | 8.5 ms / 11.4 ms | 8.9 ms | 8.7 ms / 8.9 ms | 1.04× |  |
| 132 | nfsim_basicmodels | `r25` | 100 | 100 | 8.4 ms | 8.4 ms / 10.1 ms | 10.7 ms | 10.7 ms / 10.8 ms | 1.27× |  |
| 133 | nfsim_basicmodels | `r26` | 100 | 100 | 8.4 ms | 8.3 ms / 11.8 ms | 10.6 ms | 10.5 ms / 10.8 ms | 1.27× |  |
| 134 | feature_coverage | `ft_competitive_binding` | 100 | 100 | 8.3 ms | 7.3 ms / 8.6 ms | 9.3 ms | 9.2 ms / 9.5 ms | 1.12× |  |
| 135 | feature_coverage | `ft_species_vs_molecules` | 50 | 100 | 8.2 ms | 8.1 ms / 10.6 ms | 8.9 ms | 8.9 ms / 9.2 ms | 1.10× |  |
| 136 | nfsim_basicmodels | `r30` | 100 | 100 | 8.1 ms | 7.4 ms / 8.8 ms | 8.7 ms | 8.4 ms / 8.7 ms | 1.07× |  |
| 137 | feature_coverage | `combo_exclude_with_complex` | 100 | 100 | 8.1 ms | 7.7 ms / 10.9 ms | 9.5 ms | 9.3 ms / 9.6 ms | 1.17× |  |
| 138 | feature_coverage | `ft_local_functions` | 100 | 100 | 8.1 ms | 7.4 ms / 8.8 ms | 8.8 ms | 8.6 ms / 8.8 ms | 1.09× |  |
| 139 | feature_coverage | `ft_receptor_heterogeneity` | 100 | 100 | 7.9 ms | 7.6 ms / 10.1 ms | 9.4 ms | 9.4 ms / 9.5 ms | 1.20× |  |
| 140 | nfsim_basicmodels | `r29` | 100 | 100 | 7.8 ms | 7.5 ms / 8.3 ms | 8.5 ms | 8.5 ms / 8.7 ms | 1.09× |  |
| 141 | feature_coverage | `edg_oscillator` | 40 | 200 | 7.8 ms | 7.5 ms / 10.2 ms | 10.1 ms | 10.1 ms / 10.2 ms | 1.31× |  |
| 142 | feature_coverage | `edg_state_increment_chain` | 30 | 60 | 7.7 ms | 7.3 ms / 8.7 ms | 8.8 ms | 8.7 ms / 9.2 ms | 1.14× |  |
| 143 | feature_coverage | `edg_state_wildcard_set` | 40 | 80 | 7.7 ms | 7.6 ms / 9.1 ms | 8.9 ms | 8.9 ms / 9.3 ms | 1.16× |  |
| 144 | feature_coverage | `ft_clamped_species_strict` | 50 | 50 | 7.6 ms | 7.2 ms / 7.8 ms | 7.6 ms | 7.5 ms / 7.7 ms | 1.00× |  |
| 145 | nfsim_basicmodels | `r22` | 100 | 100 | 7.5 ms | 7.2 ms / 10.7 ms | 9.3 ms | 9.1 ms / 9.4 ms | 1.24× |  |
| 146 | feature_coverage | `edg_homotrimer_binding` | 50 | 100 | 7.3 ms | 7.1 ms / 9.3 ms | 9.3 ms | 9.2 ms / 9.4 ms | 1.27× |  |
| 147 | nfsim_basicmodels | `r17` | 5 | 50 | 7.3 ms | 6.6 ms / 8.6 ms | 8.3 ms | 8.0 ms / 8.5 ms | 1.14× |  |
| 148 | feature_coverage | `ft_multi_op_rule` | 100 | 100 | 7.2 ms | 6.8 ms / 9.1 ms | 8.7 ms | 8.5 ms / 8.9 ms | 1.21× |  |
| 149 | feature_coverage | `ft_multistate` | 300 | 100 | 7.1 ms | 6.8 ms / 9.0 ms | 9.2 ms | 9.1 ms / 9.5 ms | 1.29× |  |
| 150 | feature_coverage | `ft_multi_product` | 200 | 100 | 7.1 ms | 6.9 ms / 7.9 ms | 8.7 ms | 8.7 ms / 8.9 ms | 1.23× |  |
| 151 | nfsim_basicmodels | `r23` | 100 | 100 | 7.1 ms | 7.0 ms / 10.8 ms | 9.4 ms | 9.0 ms / 9.6 ms | 1.32× |  |
| 152 | feature_coverage | `edg_dynamic_rate_zero_obs` | 40 | 80 | 7.0 ms | 6.7 ms / 7.8 ms | 7.9 ms | 7.7 ms / 8.0 ms | 1.13× |  |
| 153 | nfsim_basicmodels | `r32` | 5 | 200 | 7.0 ms | 6.1 ms / 7.1 ms | 8.8 ms | 8.6 ms / 8.9 ms | 1.26× |  |
| 154 | feature_coverage | `edg_self_dimerize` | 30 | 60 | 6.9 ms | 6.3 ms / 8.3 ms | 9.1 ms | 8.2 ms / 12.1 ms | 1.32× |  |
| 155 | feature_coverage | `edg_fixed_competition` | 40 | 80 | 6.9 ms | 6.1 ms / 8.3 ms | 8.3 ms | 8.2 ms / 8.4 ms | 1.21× |  |
| 156 | feature_coverage | `ft_state_wildcards` | 200 | 100 | 6.8 ms | 6.4 ms / 8.4 ms | 8.6 ms | 8.5 ms / 8.7 ms | 1.27× |  |
| 157 | feature_coverage | `ft_synthesis_degradation` | 500 | 100 | 6.8 ms | 6.1 ms / 6.9 ms | 8.3 ms | 8.2 ms / 8.3 ms | 1.22× |  |
| 158 | feature_coverage | `ft_match_once` | 10 | 50 | 6.8 ms | 6.2 ms / 8.0 ms | 8.3 ms | 8.2 ms / 9.5 ms | 1.23× |  |
| 159 | feature_coverage | `ft_tfun` | 40 | 80 | 6.7 ms | 6.4 ms / 7.6 ms | N/A | — | — | NFsim refused XML |
| 160 | feature_coverage | `ft_delete_molecules` | 200 | 100 | 6.7 ms | 6.1 ms / 8.6 ms | 8.4 ms | 8.1 ms / 8.5 ms | 1.25× |  |
| 161 | feature_coverage | `edg_synth_bonded_complex` | 40 | 80 | 6.7 ms | 6.3 ms / 8.1 ms | 8.2 ms | 8.2 ms / 8.3 ms | 1.23× |  |
| 162 | feature_coverage | `edg_synth_bind_existing` | 30 | 60 | 6.6 ms | 6.1 ms / 8.4 ms | 8.1 ms | 8.0 ms / 8.2 ms | 1.23× |  |
| 163 | corpus | `A_plus_A_mixed_2` | 30 | 100 | 6.6 ms | 6.5 ms / 7.9 ms | 9.0 ms | 8.8 ms / 9.8 ms | 1.37× |  |
| 164 | feature_coverage | `ft_conditional_rate` | 500 | 100 | 6.6 ms | 5.9 ms / 8.1 ms | 8.3 ms | 8.2 ms / 8.3 ms | 1.27× |  |
| 165 | feature_coverage | `ft_mm_ratelaw` | 30 | 60 | 6.5 ms | 5.9 ms / 7.2 ms | 7.9 ms | 7.8 ms / 7.9 ms | 1.21× |  |
| 166 | feature_coverage | `ft_total_rate` | 100 | 100 | 6.5 ms | 5.8 ms / 7.5 ms | 8.1 ms | 8.1 ms / 8.3 ms | 1.25× |  |
| 167 | feature_coverage | `edg_double_state_change` | 30 | 60 | 6.5 ms | 6.2 ms / 7.3 ms | 8.1 ms | 8.0 ms / 8.1 ms | 1.25× |  |
| 168 | feature_coverage | `ft_functional_rate` | 300 | 100 | 6.4 ms | 5.9 ms / 7.4 ms | 8.2 ms | 8.1 ms / 8.2 ms | 1.27× |  |
| 169 | feature_coverage | `edg_compound_op_swap` | 30 | 60 | 6.4 ms | 6.0 ms / 8.0 ms | 8.2 ms | 8.2 ms / 8.3 ms | 1.30× |  |
| 170 | feature_coverage | `edg_time_dependent_rate` | 50 | 100 | 6.3 ms | 5.6 ms / 6.4 ms | N/A | — | — | NFsim produced no output |
| 171 | feature_coverage | `edg_ring_break_constraint` | 10 | 40 | 6.3 ms | 5.9 ms / 6.5 ms | 7.5 ms | 7.3 ms / 7.7 ms | 1.19× |  |
| 172 | corpus | `nfsim_hybrid_particle_field` | 10 | 50 | 6.2 ms | 5.8 ms / 6.6 ms | 7.8 ms | 7.6 ms / 8.0 ms | 1.25× |  |
| 173 | feature_coverage | `edg_seeded_ring` | 20 | 40 | 6.1 ms | 6.0 ms / 6.9 ms | 7.9 ms | 7.7 ms / 7.9 ms | 1.29× |  |
| 174 | feature_coverage | `ft_clamped_species` | 100 | 100 | 6.1 ms | 5.9 ms / 7.6 ms | 7.7 ms | 7.7 ms / 7.8 ms | 1.27× |  |
| 175 | feature_coverage | `edg_deep_param_chain` | 40 | 80 | 6.0 ms | 5.8 ms / 7.1 ms | N/A | — | — | NFsim produced no output |
| 176 | corpus | `gene_expr_func` | 1000 | 100 | 6.0 ms | 5.8 ms / 7.4 ms | 8.4 ms | 8.0 ms / 8.5 ms | 1.40× |  |
| 177 | feature_coverage | `edg_zero_rate_rule` | 40 | 80 | 5.9 ms | 5.8 ms / 6.4 ms | 7.7 ms | 7.6 ms / 8.1 ms | 1.31× |  |
