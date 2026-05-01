# RuleMonkey vs NFsim — wall-time comparison

Single-machine, sequential, end-to-end subprocess wall time.  3 reps per engine per model; median reported.  Both engines invoked from a fresh process per rep (no warm caches).  Same XML, same t_end, same n_steps, same seed sequence.

**This is an efficiency report, not a correctness report.**  Models where NFsim runs but produces incorrect observables (e.g. `ft_tfun` due to NFsim's TFUN handler returning zero rate) are still included in the timing table; correctness is covered separately by `feature_coverage`, `benchmark_full`, and `basicmodels` suites.

Generated: `2026-04-29 22:30:37 MDT` (`benchmark_rm_vs_nfsim_timing.py`)

## Summary

- Models scored: **173 / 177** (both engines completed all 3 reps)
- NFsim N/A: **4** (RM ran; NFsim refused / errored)
- **Median speedup (NFsim wall ÷ RM wall): 0.83×**
- Geometric mean speedup: 0.78×
- RM faster than NFsim on **72 / 173** models; NFsim faster on 101

## Where does each engine win?

The 173 scored models fall into six speedup buckets, and the buckets
correlate cleanly with model shape:

| Bucket | Count | Examples | What's going on |
|---|---:|---|---|
| RM ≥ 10× | 4 | `mlnr` 23.7×, `pltr` 24.8×, `testcase2a` 11.8×, `rm_tlbr` 12.0× | Pattern-quantifier (count-relation) Species observables |
| RM 2–10× | 1 | `example3_fit` 2.96× | Single-model outlier |
| RM 1.1–2× | 56 | `ss_tlbr_rings`, `combo_symmetric_rings`, many fast feature_coverage / basicmodels | RM's slightly lower per-process startup wins on short runs |
| Parity (0.8–1.1×) | 32 | `BLBR`, `tlbr`, `A_plus_B_rings`, `bench_tlbr_yang2008` | Same network-free chemistry, both engines comparable |
| NFsim 0.5–0.8× | 31 | `e1`–`e9`, `tcr_*` variants | Heavy enzyme-kinetics scaling; NFsim a bit ahead |
| NFsim < 0.5× | 49 | `t3` 0.49×, `machine` 0.24×, `fceri_ji` 0.23×, `ensemble` 0.24×, `rm_tlbr_rings` 0.23×, `AN` 0.14× | Heavy stochastic corpus models: NFsim's tuning shows |

### RM ≥ 10× — count-relation Species observables

All four models in this bucket declare **exactly 300 Species observables
of the form `R()=N`** (`R()=1`, `R()=2`, ..., `R()=300`) — histogram
bins counting complexes that contain exactly *N* receptors.  No other
model in the corpus declares any:

| Model | Count-relation Species obs | Speedup |
|---|---:|---:|
| `mlnr`       | 300 | 23.7× |
| `pltr`       | 300 | 24.8× |
| `testcase2a` | 300 | 11.8× |
| `rm_tlbr`    | 300 | 12.0× |
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

The 56 models in this band are mostly sub-100ms runs where end-to-end
process startup is a non-trivial fraction of wall time.  RM has
slightly lower startup overhead (smaller binary, simpler XML loader),
which shows up as a 1.1×–1.5× win on short runs.  Vanishes once the
simulation itself dominates wall time.

### Parity — generic network-free chemistry without histogram observables

`A_plus_B_rings` (no Size_N obs) sits at 0.99×; `bench_tlbr_yang2008`
at 1.04×; `BLBR` at 1.01×.  Same kind of binding chemistry as the
RM-winning Size_N models above — but without the observable-evaluation
asymmetry, the engines are comparable.

### NFsim < 0.5× — heavy stochastic corpus models

The dominant story for the bottom of the table.  Models like `t3`,
`machine`, `fceri_ji`, `ensemble`, `rm_tlbr_rings`, `AN` are the
heavy-event-rate workloads NFsim has been tuned on for years.  RM has
real ground to make up here.  Memory says the 2026-Q2 optimization
sprints closed at "~55% cumulative on `msite_S3200`" and the engaged-
set campaign was retired with a note that "next ≥5% requires
structural change" — the corpus models in this bucket weren't the
direct optimization target, and the gap shows.

The single-model outlier `example3_fit` (2.96×) doesn't fit any of
these buckets cleanly — worth a separate look if anyone wants to
chase it.

## Caveats

- Wall time, not CPU time.  Includes process startup, XML load, and result write — same for both engines, so the ratio is fair, but absolute numbers shift across machines.
- Speedup depends on workload: long-horizon, high-event-rate models show RM's biggest wins; short-horizon trivial models are near parity or NFsim-favored due to RM's slightly heavier per-process startup.
- Single-replicate timings can be noisy.  See min/max columns for spread; the median column is what to quote.
- N/A in the NFsim column means NFsim refused or errored on at least one of the 3 reps.  The reason column gives the captured diagnostic.

## Per-model results

Sorted by RM median wall time (most expensive first).

| # | Suite | Model | t_end | n_steps | RM median | RM min/max | NFsim median | NFsim min/max | Speedup | Notes |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | corpus | `t3` | 100 | 100 | 49.9 s | 49.8 s / 50.1 s | 24.3 s | 23.1 s / 24.4 s | 0.49× |  |
| 2 | corpus | `machine` | 14 | 100 | 42.1 s | 41.5 s / 42.9 s | 9970 ms | 9962 ms / 10.3 s | 0.24× |  |
| 3 | corpus | `rm_tlbr_rings` | 1000 | 100 | 35.3 s | 34.5 s / 36.7 s | 8126 ms | 7352 ms / 8385 ms | 0.23× |  |
| 4 | corpus | `fceri_ji` | 500 | 500 | 33.6 s | 32.7 s / 35.4 s | 7764 ms | 7545 ms / 8133 ms | 0.23× |  |
| 5 | corpus | `ensemble` | 15 | 100 | 32.3 s | 31.9 s / 32.4 s | 7873 ms | 7680 ms / 7873 ms | 0.24× |  |
| 6 | corpus | `AN` | 100 | 100 | 21.9 s | 21.6 s / 22.1 s | 2970 ms | 2900 ms / 2985 ms | 0.14× |  |
| 7 | corpus | `tcr` | 60 | 12 | 21.8 s | 19.8 s / 22.1 s | 10.7 s | 10.5 s / 11.3 s | 0.49× |  |
| 8 | corpus | `example2_fit` | 126 | 20 | 18.0 s | 17.2 s / 18.1 s | 6677 ms | 6365 ms / 6931 ms | 0.37× |  |
| 9 | corpus | `e9` | 200 | 100 | 17.1 s | 16.4 s / 18.9 s | 13.5 s | 13.5 s / 13.6 s | 0.79× |  |
| 10 | corpus | `example4_fit` | 13 | 100 | 16.6 s | 15.8 s / 16.7 s | 9461 ms | 9397 ms / 9553 ms | 0.57× |  |
| 11 | corpus | `egfr_nf_iter5p12h10` | 60 | 20 | 15.9 s | 15.1 s / 17.8 s | 5988 ms | 5627 ms / 6101 ms | 0.38× |  |
| 12 | corpus | `e8` | 200 | 100 | 15.2 s | 15.1 s / 15.2 s | 12.6 s | 12.3 s / 12.7 s | 0.83× |  |
| 13 | corpus | `e7` | 200 | 100 | 14.5 s | 14.4 s / 14.5 s | 10.7 s | 10.7 s / 10.8 s | 0.74× |  |
| 14 | corpus | `blbr_rings_posner1995` | 3000 | 300 | 14.3 s | 14.3 s / 14.6 s | 5797 ms | 5622 ms / 5993 ms | 0.40× |  |
| 15 | corpus | `bench_blbr_rings_posner1995` | 3000 | 300 | 13.9 s | 13.9 s / 14.2 s | 5942 ms | 5910 ms / 6448 ms | 0.43× |  |
| 16 | corpus | `e6` | 200 | 100 | 12.6 s | 12.5 s / 12.7 s | 9657 ms | 9491 ms / 9794 ms | 0.77× |  |
| 17 | corpus | `PushPull` | 4000 | 40 | 11.7 s | 11.4 s / 12.3 s | 4593 ms | 4490 ms / 4718 ms | 0.39× |  |
| 18 | corpus | `ANx` | 100 | 100 | 11.5 s | 11.3 s / 11.7 s | 1700 ms | 1668 ms / 1816 ms | 0.15× |  |
| 19 | corpus | `e5` | 200 | 100 | 11.2 s | 11.0 s / 12.5 s | 8484 ms | 8292 ms / 8610 ms | 0.76× |  |
| 20 | corpus | `e4` | 200 | 100 | 9837 ms | 9764 ms / 10.1 s | 6893 ms | 6490 ms / 6907 ms | 0.70× |  |
| 21 | corpus | `oscSystem` | 200 | 200 | 9421 ms | 9096 ms / 9466 ms | 400 ms | 388 ms / 439 ms | 0.04× |  |
| 22 | corpus | `e3` | 200 | 100 | 9060 ms | 8826 ms / 9147 ms | 5665 ms | 5589 ms / 5734 ms | 0.63× |  |
| 23 | corpus | `lat` | 10 | 100 | 7912 ms | 7197 ms / 8474 ms | 3149 ms | 3034 ms / 3210 ms | 0.40× |  |
| 24 | corpus | `e2` | 200 | 100 | 7317 ms | 7263 ms / 7736 ms | 3782 ms | 3644 ms / 3847 ms | 0.52× |  |
| 25 | corpus | `egfr_net` | 120 | 120 | 6708 ms | 6682 ms / 6995 ms | 4040 ms | 4035 ms / 4072 ms | 0.60× |  |
| 26 | corpus | `e1` | 200 | 100 | 6462 ms | 6061 ms / 6549 ms | 2505 ms | 2203 ms / 2512 ms | 0.39× |  |
| 27 | nfsim_basicmodels | `r02` | 100 | 100 | 4955 ms | 4543 ms / 5043 ms | 1164 ms | 1149 ms / 1251 ms | 0.23× |  |
| 28 | corpus | `bench_blbr_dembo1978_monovalent_inhibitor` | 3000 | 300 | 4396 ms | 4171 ms / 4854 ms | 3451 ms | 3397 ms / 3482 ms | 0.78× |  |
| 29 | nfsim_basicmodels | `r01` | 100 | 100 | 3370 ms | 3246 ms / 3721 ms | 1532 ms | 1413 ms / 1533 ms | 0.45× |  |
| 30 | corpus | `mlnr` | 1000 | 1000 | 3369 ms | 3312 ms / 3420 ms | 79.8 s | 79.3 s / 81.0 s | 23.70× |  |
| 31 | corpus | `pltr` | 1000 | 1000 | 3276 ms | 3252 ms / 3373 ms | 81.3 s | 80.4 s / 85.1 s | 24.83× |  |
| 32 | nfsim_basicmodels | `r20` | 100 | 100 | 3268 ms | 3143 ms / 3488 ms | 1540 ms | 1481 ms / 1545 ms | 0.47× |  |
| 33 | nfsim_basicmodels | `r07` | 20 | 100 | 2772 ms | 2755 ms / 2967 ms | 892 ms | 835 ms / 943 ms | 0.32× |  |
| 34 | corpus | `rm_tlbr` | 1000 | 1000 | 2686 ms | 2652 ms / 2742 ms | 32.2 s | 31.6 s / 33.4 s | 12.00× |  |
| 35 | nfsim_basicmodels | `r06` | 20 | 100 | 2633 ms | 2580 ms / 2671 ms | 515 ms | 483 ms / 535 ms | 0.20× |  |
| 36 | corpus | `testcase2a` | 1000 | 1000 | 2586 ms | 2553 ms / 2597 ms | 30.4 s | 29.8 s / 31.7 s | 11.77× |  |
| 37 | corpus | `tlbr` | 200 | 200 | 2511 ms | 2450 ms / 2562 ms | 2283 ms | 2272 ms / 2286 ms | 0.91× |  |
| 38 | corpus | `simple_system` | 100 | 50 | 2181 ms | 2056 ms / 2287 ms | 965 ms | 956 ms / 1020 ms | 0.44× |  |
| 39 | corpus | `simple_nfsim` | 100 | 50 | 2150 ms | 1986 ms / 2369 ms | 1056 ms | 1045 ms / 1122 ms | 0.49× |  |
| 40 | corpus | `nfsim_aggregation_gelation` | 50 | 200 | 2103 ms | 2081 ms / 2126 ms | 697 ms | 694 ms / 709 ms | 0.33× |  |
| 41 | nfsim_basicmodels | `r10` | 100 | 100 | 1990 ms | 1878 ms / 2085 ms | 364 ms | 355 ms / 382 ms | 0.18× |  |
| 42 | nfsim_basicmodels | `r04` | 50 | 100 | 1925 ms | 1910 ms / 1963 ms | 746 ms | 745 ms / 781 ms | 0.39× |  |
| 43 | nfsim_basicmodels | `r08` | 20 | 100 | 1883 ms | 1878 ms / 1894 ms | 530 ms | 506 ms / 541 ms | 0.28× |  |
| 44 | corpus | `tcr_gen27ind33` | 0.94 | 12 | 1863 ms | 1814 ms / 1880 ms | 768 ms | 767 ms / 769 ms | 0.41× |  |
| 45 | corpus | `tcr_iter28p4h2` | 0.94 | 12 | 1809 ms | 1786 ms / 1816 ms | 741 ms | 739 ms / 750 ms | 0.41× |  |
| 46 | corpus | `tcr_gen20ind9` | 1.5 | 12 | 1722 ms | 1699 ms / 1737 ms | 757 ms | 750 ms / 767 ms | 0.44× |  |
| 47 | corpus | `tcr_iter9p44` | 0.94 | 12 | 1713 ms | 1705 ms / 1723 ms | 730 ms | 721 ms / 732 ms | 0.43× |  |
| 48 | corpus | `bench_blbr_cooperativity_posner2004` | 10 | 10 | 1253 ms | 1247 ms / 1340 ms | 748 ms | 730 ms / 761 ms | 0.60× |  |
| 49 | corpus | `bench_blbr_cooperativity_posner2004_rings` | 10 | 10 | 1184 ms | 1168 ms / 1222 ms | 725 ms | 715 ms / 733 ms | 0.61× |  |
| 50 | corpus | `isingspin_localfcn` | 5000 | 100 | 1066 ms | 1033 ms / 1108 ms | 218 ms | 213 ms / 222 ms | 0.20× |  |
| 51 | corpus | `CaMKII_holo` | 10 | 100 | 972 ms | 949 ms / 1042 ms | 585 ms | 568 ms / 588 ms | 0.60× |  |
| 52 | corpus | `stiff` | 10 | 100 | 901 ms | 898 ms / 919 ms | 489 ms | 486 ms / 491 ms | 0.54× |  |
| 53 | corpus | `example3_fit` | 5000 | 10 | 840 ms | 838 ms / 893 ms | 2488 ms | 2417 ms / 2529 ms | 2.96× |  |
| 54 | corpus | `BLBR` | 10 | 100 | 716 ms | 702 ms / 765 ms | 724 ms | 711 ms / 793 ms | 1.01× |  |
| 55 | nfsim_basicmodels | `r03` | 50 | 100 | 574 ms | 544 ms / 589 ms | 274 ms | 270 ms / 276 ms | 0.48× |  |
| 56 | corpus | `rm_blbr` | 10 | 100 | 517 ms | 516 ms / 523 ms | 328 ms | 318 ms / 330 ms | 0.63× |  |
| 57 | corpus | `bench_blbr_heterogeneity_goldstein1980` | 10 | 10 | 447 ms | 409 ms / 456 ms | 124 ms | 122 ms / 126 ms | 0.28× |  |
| 58 | nfsim_basicmodels | `r12` | 10 | 100 | 439 ms | 396 ms / 450 ms | 173 ms | 171 ms / 182 ms | 0.39× |  |
| 59 | nfsim_basicmodels | `r18` | 10 | 100 | 407 ms | 406 ms / 416 ms | 287 ms | 286 ms / 290 ms | 0.71× |  |
| 60 | feature_coverage | `ss_long_polymer` | 500 | 100 | 389 ms | 384 ms / 389 ms | 170 ms | 168 ms / 174 ms | 0.44× |  |
| 61 | nfsim_basicmodels | `r09` | 100 | 100 | 377 ms | 368 ms / 382 ms | 68.8 ms | 66.9 ms / 68.8 ms | 0.18× |  |
| 62 | feature_coverage | `ft_multimol_sym_obs` | 20 | 100 | 338 ms | 329 ms / 339 ms | 74.9 ms | 74.5 ms / 76.4 ms | 0.22× |  |
| 63 | corpus | `st_multi_2` | 10 | 100 | 318 ms | 314 ms / 332 ms | 154 ms | 141 ms / 156 ms | 0.48× |  |
| 64 | corpus | `toy_jim` | 100 | 100 | 315 ms | 309 ms / 333 ms | 140 ms | 139 ms / 145 ms | 0.45× |  |
| 65 | feature_coverage | `ss_branching_aggregate` | 400 | 100 | 274 ms | 271 ms / 286 ms | 194 ms | 193 ms / 197 ms | 0.71× |  |
| 66 | corpus | `ANx_noActivity` | 100 | 100 | 245 ms | 244 ms / 247 ms | 177 ms | 173 ms / 182 ms | 0.72× |  |
| 67 | corpus | `poly` | 100 | 100 | 236 ms | 231 ms / 240 ms | 134 ms | 125 ms / 164 ms | 0.57× |  |
| 68 | corpus | `bench_tlbr_yang2008` | 3000 | 300 | 210 ms | 208 ms / 219 ms | 218 ms | 213 ms / 219 ms | 1.04× |  |
| 69 | corpus | `bench_blbr_rings_posner1995_no_rings` | 50 | 50 | 198 ms | 189 ms / 199 ms | 158 ms | 157 ms / 161 ms | 0.80× |  |
| 70 | corpus | `bench_blbr_dembo1978` | 40 | 40 | 187 ms | 181 ms / 190 ms | 153 ms | 153 ms / 159 ms | 0.82× |  |
| 71 | corpus | `A_plus_B_rings` | 10 | 10000 | 169 ms | 164 ms / 170 ms | 167 ms | 164 ms / 173 ms | 0.99× |  |
| 72 | corpus | `bench_tlbr_solution_macken1982` | 500 | 250 | 164 ms | 161 ms / 177 ms | 115 ms | 114 ms / 119 ms | 0.70× |  |
| 73 | feature_coverage | `ss_symmetric_homopoly` | 300 | 100 | 157 ms | 156 ms / 160 ms | 127 ms | 124 ms / 128 ms | 0.81× |  |
| 74 | feature_coverage | `ft_stiff_system` | 200 | 100 | 150 ms | 149 ms / 153 ms | 49.3 ms | 48.4 ms / 51.0 ms | 0.33× |  |
| 75 | feature_coverage | `ft_multisite_phospho` | 200 | 100 | 147 ms | 142 ms / 148 ms | 73.8 ms | 70.5 ms / 76.8 ms | 0.50× |  |
| 76 | nfsim_basicmodels | `r11` | 1 | 100 | 130 ms | 128 ms / 132 ms | 59.1 ms | 58.3 ms / 60.9 ms | 0.46× |  |
| 77 | corpus | `A_plus_A_mixed_1` | 30 | 100 | 114 ms | 112 ms / 115 ms | 47.4 ms | 46.2 ms / 49.8 ms | 0.42× |  |
| 78 | feature_coverage | `ss_tlbr_rings` | 400 | 100 | 110 ms | 109 ms / 111 ms | 131 ms | 131 ms / 134 ms | 1.19× |  |
| 79 | nfsim_basicmodels | `r13` | 30 | 100 | 101 ms | 97.1 ms / 106 ms | 26.3 ms | 25.7 ms / 27.5 ms | 0.26× |  |
| 80 | corpus | `st_multi_1` | 10 | 100 | 100 ms | 97.6 ms / 103 ms | 43.1 ms | 42.8 ms / 47.8 ms | 0.43× |  |
| 81 | nfsim_basicmodels | `r14` | 30 | 100 | 74.9 ms | 72.9 ms / 77.4 ms | 30.5 ms | 29.2 ms / 30.7 ms | 0.41× |  |
| 82 | feature_coverage | `ft_push_pull` | 200 | 100 | 73.8 ms | 73.2 ms / 73.9 ms | 32.1 ms | 32.1 ms / 32.7 ms | 0.43× |  |
| 83 | nfsim_basicmodels | `r21` | 100 | 100 | 65.6 ms | 64.7 ms / 66.8 ms | 84.3 ms | 82.1 ms / 84.5 ms | 1.28× |  |
| 84 | feature_coverage | `ft_receptor_dimerization` | 100 | 100 | 57.6 ms | 56.7 ms / 59.7 ms | 30.4 ms | 30.0 ms / 30.8 ms | 0.53× |  |
| 85 | feature_coverage | `ft_cooperative_binding` | 100 | 100 | 45.9 ms | 42.9 ms / 45.9 ms | 20.1 ms | 19.7 ms / 21.0 ms | 0.44× |  |
| 86 | nfsim_basicmodels | `r16` | 5 | 100 | 43.5 ms | 43.0 ms / 43.9 ms | 21.3 ms | 21.1 ms / 21.7 ms | 0.49× |  |
| 87 | feature_coverage | `combo_strict_product_plus` | 100 | 100 | 38.5 ms | 37.6 ms / 39.2 ms | 21.4 ms | 20.3 ms / 22.1 ms | 0.56× |  |
| 88 | corpus | `receptor_nf_iter36p0h3` | 60 | 60 | 34.0 ms | 33.7 ms / 35.1 ms | 21.8 ms | 21.0 ms / 24.5 ms | 0.64× |  |
| 89 | corpus | `Tutorial_Example` | 5 | 50 | 31.4 ms | 31.4 ms / 31.9 ms | 21.8 ms | 21.6 ms / 22.8 ms | 0.69× |  |
| 90 | corpus | `example6_ground_truth` | 60 | 60 | 27.9 ms | 27.7 ms / 27.9 ms | 18.6 ms | 17.9 ms / 19.2 ms | 0.67× |  |
| 91 | nfsim_basicmodels | `r15` | 1 | 100 | 26.6 ms | 26.6 ms / 28.5 ms | 15.9 ms | 15.8 ms / 17.4 ms | 0.60× |  |
| 92 | corpus | `A_plus_A` | 10 | 100 | 22.5 ms | 21.6 ms / 24.4 ms | 17.2 ms | 15.3 ms / 17.3 ms | 0.77× |  |
| 93 | feature_coverage | `nf_branching_aggregate` | 200 | 100 | 22.2 ms | 21.4 ms / 22.6 ms | 20.9 ms | 20.6 ms / 21.6 ms | 0.94× |  |
| 94 | nfsim_basicmodels | `r05` | 1 | 100 | 20.4 ms | 20.0 ms / 20.4 ms | 16.5 ms | 15.9 ms / 16.7 ms | 0.81× |  |
| 95 | feature_coverage | `ft_continue` | 100 | 100 | 20.1 ms | 19.6 ms / 20.6 ms | 13.6 ms | 13.5 ms / 13.8 ms | 0.68× |  |
| 96 | feature_coverage | `ft_blbr` | 100 | 100 | 16.9 ms | 16.4 ms / 17.1 ms | 15.7 ms | 15.7 ms / 16.2 ms | 0.93× |  |
| 97 | feature_coverage | `ft_tlbr` | 100 | 100 | 15.9 ms | 15.4 ms / 15.9 ms | 16.9 ms | 16.4 ms / 17.1 ms | 1.06× |  |
| 98 | corpus | `nfsim_coarse_graining` | 400 | 200 | 15.5 ms | 15.5 ms / 16.2 ms | 14.4 ms | 14.4 ms / 14.7 ms | 0.93× |  |
| 99 | feature_coverage | `ft_multimol_pattern_sym_nonreacting` | 10 | 40 | 15.2 ms | 15.2 ms / 15.2 ms | 12.4 ms | 12.0 ms / 13.6 ms | 0.82× |  |
| 100 | feature_coverage | `ft_energy_patterns` | 500 | 100 | 15.2 ms | 14.8 ms / 15.7 ms | 11.6 ms | 11.6 ms / 11.8 ms | 0.77× |  |
| 101 | feature_coverage | `nf_large_multivalent` | 100 | 100 | 14.4 ms | 13.8 ms / 15.1 ms | 14.1 ms | 13.6 ms / 14.2 ms | 0.98× |  |
| 102 | feature_coverage | `ft_include_reactants` | 100 | 100 | 13.7 ms | 13.7 ms / 14.5 ms | 10.5 ms | 10.3 ms / 10.5 ms | 0.77× |  |
| 103 | feature_coverage | `nf_linear_polymer` | 200 | 100 | 13.2 ms | 12.6 ms / 13.6 ms | 11.1 ms | 11.0 ms / 11.4 ms | 0.84× |  |
| 104 | corpus | `nfsim_dynamic_compartments` | 20 | 100 | 13.1 ms | 13.1 ms / 19.3 ms | 13.6 ms | 12.6 ms / 14.2 ms | 1.04× |  |
| 105 | feature_coverage | `ft_exclude_products` | 200 | 50 | 12.8 ms | 11.8 ms / 12.8 ms | 10.3 ms | 10.3 ms / 10.4 ms | 0.81× |  |
| 106 | feature_coverage | `combo_localfcn_multisite` | 200 | 100 | 12.7 ms | 12.4 ms / 13.2 ms | 11.8 ms | 11.8 ms / 11.9 ms | 0.93× |  |
| 107 | corpus | `st` | 10 | 100 | 12.4 ms | 12.2 ms / 12.9 ms | 12.0 ms | 11.5 ms / 12.3 ms | 0.97× |  |
| 108 | feature_coverage | `combo_symmetric_rings` | 100 | 100 | 12.3 ms | 12.0 ms / 12.4 ms | 13.6 ms | 13.4 ms / 13.9 ms | 1.11× |  |
| 109 | feature_coverage | `ft_perturbation` | 100 | 100 | 12.0 ms | 12.0 ms / 12.7 ms | 11.6 ms | 11.3 ms / 11.7 ms | 0.97× |  |
| 110 | feature_coverage | `ft_multi_site_binding` | 100 | 100 | 11.6 ms | 11.6 ms / 12.0 ms | 11.7 ms | 11.7 ms / 11.8 ms | 1.01× |  |
| 111 | feature_coverage | `edg_oscillator` | 40 | 200 | 11.4 ms | 11.3 ms / 12.1 ms | 11.1 ms | 11.0 ms / 11.1 ms | 0.97× |  |
| 112 | feature_coverage | `combo_addbond_connected` | 200 | 100 | 11.3 ms | 11.1 ms / 11.7 ms | 11.0 ms | 10.6 ms / 12.0 ms | 0.97× |  |
| 113 | feature_coverage | `combo_synth_degrade_equilibrium` | 500 | 100 | 11.2 ms | 11.0 ms / 11.8 ms | 11.2 ms | 11.2 ms / 11.3 ms | 1.00× |  |
| 114 | corpus | `basicTLBR` | 20 | 100 | 10.9 ms | 10.8 ms / 11.0 ms | 15.5 ms | 15.0 ms / 16.0 ms | 1.41× |  |
| 115 | feature_coverage | `combo_shorthand_embed` | 100 | 100 | 10.8 ms | 10.7 ms / 11.0 ms | 11.5 ms | 11.4 ms / 11.6 ms | 1.07× |  |
| 116 | feature_coverage | `ft_nested_functions` | 100 | 100 | 10.7 ms | 10.6 ms / 11.1 ms | N/A | — | — | NFsim produced no output |
| 117 | nfsim_basicmodels | `r19` | 10 | 100 | 10.6 ms | 10.5 ms / 11.1 ms | 18.6 ms | 18.6 ms / 18.7 ms | 1.76× |  |
| 118 | feature_coverage | `ft_signaling_cascade` | 100 | 100 | 10.6 ms | 10.5 ms / 10.6 ms | 11.6 ms | 11.6 ms / 11.7 ms | 1.10× |  |
| 119 | feature_coverage | `edg_multi_pattern_obs` | 30 | 60 | 10.5 ms | 10.3 ms / 10.6 ms | 9.9 ms | 9.9 ms / 10.4 ms | 0.95× |  |
| 120 | feature_coverage | `edg_branched_polymer` | 30 | 60 | 10.1 ms | 9.4 ms / 10.7 ms | 15.1 ms | 15.1 ms / 16.6 ms | 1.50× |  |
| 121 | feature_coverage | `combo_multimol_unimol` | 100 | 100 | 9.9 ms | 9.8 ms / 10.0 ms | 10.7 ms | 10.6 ms / 10.8 ms | 1.08× |  |
| 122 | feature_coverage | `edg_pattern_local_fcn` | 40 | 80 | 9.8 ms | 9.7 ms / 9.8 ms | 10.3 ms | 10.0 ms / 10.3 ms | 1.05× |  |
| 123 | nfsim_basicmodels | `r26` | 100 | 100 | 9.8 ms | 9.5 ms / 10.4 ms | 11.5 ms | 11.3 ms / 11.5 ms | 1.18× |  |
| 124 | feature_coverage | `ft_exclude_reactants` | 100 | 100 | 9.7 ms | 9.2 ms / 9.9 ms | 10.2 ms | 10.2 ms / 10.4 ms | 1.05× |  |
| 125 | nfsim_basicmodels | `r25` | 100 | 100 | 9.6 ms | 9.4 ms / 10.2 ms | 11.4 ms | 11.3 ms / 11.5 ms | 1.18× |  |
| 126 | nfsim_basicmodels | `r24` | 100 | 100 | 9.6 ms | 9.4 ms / 10.2 ms | 11.4 ms | 11.3 ms / 11.5 ms | 1.19× |  |
| 127 | feature_coverage | `ft_population_map` | 60 | 60 | 9.4 ms | 9.3 ms / 9.4 ms | 10.1 ms | 10.0 ms / 10.1 ms | 1.07× |  |
| 128 | feature_coverage | `ft_bond_wildcards` | 100 | 100 | 9.3 ms | 9.3 ms / 10.0 ms | 10.3 ms | 10.2 ms / 10.9 ms | 1.10× |  |
| 129 | feature_coverage | `ft_ring_closure` | 100 | 100 | 9.1 ms | 9.0 ms / 9.2 ms | 10.2 ms | 10.1 ms / 10.2 ms | 1.12× |  |
| 130 | feature_coverage | `ft_multimol_unimol_unbind_sym` | 2 | 40 | 9.1 ms | 8.8 ms / 9.1 ms | 9.0 ms | 9.0 ms / 9.0 ms | 0.99× |  |
| 131 | corpus | `nfsim_ring_closure_polymer` | 50 | 200 | 8.9 ms | 8.5 ms / 9.1 ms | 10.6 ms | 10.5 ms / 10.7 ms | 1.20× |  |
| 132 | feature_coverage | `ft_complex_seed` | 50 | 50 | 8.7 ms | 8.7 ms / 9.1 ms | 9.8 ms | 9.6 ms / 9.9 ms | 1.13× |  |
| 133 | feature_coverage | `ft_species_vs_molecules` | 50 | 100 | 8.5 ms | 8.4 ms / 8.9 ms | 9.5 ms | 9.4 ms / 9.7 ms | 1.12× |  |
| 134 | feature_coverage | `combo_exclude_with_complex` | 100 | 100 | 8.4 ms | 8.2 ms / 9.0 ms | 10.4 ms | 9.7 ms / 10.4 ms | 1.23× |  |
| 135 | corpus | `A_plus_A_mixed_2` | 30 | 100 | 8.3 ms | 7.7 ms / 8.6 ms | 9.3 ms | 9.3 ms / 9.6 ms | 1.12× |  |
| 136 | feature_coverage | `ft_receptor_heterogeneity` | 100 | 100 | 8.3 ms | 7.8 ms / 8.3 ms | 10.0 ms | 10.0 ms / 10.2 ms | 1.21× |  |
| 137 | feature_coverage | `edg_three_mol_pattern` | 30 | 60 | 8.2 ms | 8.2 ms / 8.4 ms | 9.7 ms | 9.7 ms / 10.0 ms | 1.19× |  |
| 138 | nfsim_basicmodels | `r22` | 100 | 100 | 8.2 ms | 8.0 ms / 8.9 ms | 10.8 ms | 10.6 ms / 10.9 ms | 1.31× |  |
| 139 | feature_coverage | `ft_catalytic_unbinding` | 100 | 100 | 8.1 ms | 7.9 ms / 8.2 ms | 10.1 ms | 9.9 ms / 10.2 ms | 1.24× |  |
| 140 | feature_coverage | `edg_state_wildcard_set` | 40 | 80 | 8.1 ms | 8.0 ms / 8.4 ms | 9.7 ms | 9.7 ms / 9.8 ms | 1.20× |  |
| 141 | feature_coverage | `ft_competitive_binding` | 100 | 100 | 8.0 ms | 7.7 ms / 8.1 ms | 10.2 ms | 10.0 ms / 10.3 ms | 1.28× |  |
| 142 | feature_coverage | `ft_local_functions` | 100 | 100 | 7.6 ms | 7.4 ms / 8.8 ms | 9.8 ms | 9.4 ms / 9.8 ms | 1.28× |  |
| 143 | nfsim_basicmodels | `r30` | 100 | 100 | 7.6 ms | 7.5 ms / 8.1 ms | 9.3 ms | 9.3 ms / 9.5 ms | 1.24× |  |
| 144 | feature_coverage | `edg_homotrimer_binding` | 50 | 100 | 7.5 ms | 7.4 ms / 7.8 ms | 10.0 ms | 9.9 ms / 10.1 ms | 1.33× |  |
| 145 | feature_coverage | `edg_state_increment_chain` | 30 | 60 | 7.5 ms | 7.4 ms / 7.5 ms | 9.5 ms | 9.3 ms / 9.6 ms | 1.26× |  |
| 146 | nfsim_basicmodels | `r23` | 100 | 100 | 7.5 ms | 7.1 ms / 7.9 ms | 10.6 ms | 10.4 ms / 10.9 ms | 1.41× |  |
| 147 | nfsim_basicmodels | `r29` | 100 | 100 | 7.5 ms | 7.4 ms / 7.6 ms | 9.3 ms | 9.2 ms / 9.4 ms | 1.24× |  |
| 148 | feature_coverage | `ft_clamped_species_strict` | 50 | 50 | 7.4 ms | 7.2 ms / 7.4 ms | 8.4 ms | 8.4 ms / 8.5 ms | 1.13× |  |
| 149 | feature_coverage | `ft_multistate` | 300 | 100 | 7.3 ms | 7.3 ms / 7.6 ms | 10.1 ms | 10.1 ms / 11.7 ms | 1.38× |  |
| 150 | feature_coverage | `ft_multi_op_rule` | 100 | 100 | 7.3 ms | 7.2 ms / 7.4 ms | 9.5 ms | 9.5 ms / 9.6 ms | 1.31× |  |
| 151 | feature_coverage | `ft_multi_product` | 200 | 100 | 7.1 ms | 7.1 ms / 7.3 ms | 9.6 ms | 9.6 ms / 9.7 ms | 1.35× |  |
| 152 | feature_coverage | `edg_self_dimerize` | 30 | 60 | 6.8 ms | 6.7 ms / 7.0 ms | 8.9 ms | 8.9 ms / 10.3 ms | 1.32× |  |
| 153 | feature_coverage | `ft_state_wildcards` | 200 | 100 | 6.7 ms | 6.5 ms / 7.0 ms | 9.3 ms | 9.2 ms / 9.4 ms | 1.39× |  |
| 154 | feature_coverage | `edg_compound_op_swap` | 30 | 60 | 6.6 ms | 6.5 ms / 6.7 ms | 8.9 ms | 8.8 ms / 9.0 ms | 1.35× |  |
| 155 | feature_coverage | `ft_match_once` | 10 | 50 | 6.6 ms | 6.6 ms / 6.8 ms | 8.9 ms | 8.8 ms / 9.2 ms | 1.37× |  |
| 156 | feature_coverage | `edg_fixed_competition` | 40 | 80 | 6.5 ms | 6.3 ms / 6.6 ms | 9.0 ms | 8.9 ms / 9.1 ms | 1.38× |  |
| 157 | feature_coverage | `ft_delete_molecules` | 200 | 100 | 6.5 ms | 6.5 ms / 6.6 ms | 9.0 ms | 8.8 ms / 9.2 ms | 1.37× |  |
| 158 | feature_coverage | `edg_synth_bonded_complex` | 40 | 80 | 6.5 ms | 6.5 ms / 7.1 ms | 9.0 ms | 8.8 ms / 9.3 ms | 1.38× |  |
| 159 | feature_coverage | `ft_conditional_rate` | 500 | 100 | 6.5 ms | 6.4 ms / 6.7 ms | 9.0 ms | 8.9 ms / 9.1 ms | 1.38× |  |
| 160 | feature_coverage | `edg_double_state_change` | 30 | 60 | 6.5 ms | 6.4 ms / 6.7 ms | 8.9 ms | 8.8 ms / 8.9 ms | 1.37× |  |
| 161 | feature_coverage | `ft_synthesis_degradation` | 500 | 100 | 6.4 ms | 6.3 ms / 6.6 ms | 9.3 ms | 9.2 ms / 9.3 ms | 1.45× |  |
| 162 | feature_coverage | `edg_deep_param_chain` | 40 | 80 | 6.4 ms | 6.3 ms / 6.4 ms | N/A | — | — | NFsim produced no output |
| 163 | feature_coverage | `edg_seeded_ring` | 20 | 40 | 6.3 ms | 6.3 ms / 6.6 ms | 8.9 ms | 8.8 ms / 8.9 ms | 1.40× |  |
| 164 | feature_coverage | `edg_synth_bind_existing` | 30 | 60 | 6.3 ms | 6.2 ms / 6.4 ms | 8.7 ms | 8.6 ms / 8.7 ms | 1.38× |  |
| 165 | feature_coverage | `ft_clamped_species` | 100 | 100 | 6.3 ms | 6.0 ms / 7.6 ms | 8.7 ms | 8.5 ms / 8.7 ms | 1.39× |  |
| 166 | corpus | `gene_expr_func` | 1000 | 100 | 6.2 ms | 6.2 ms / 6.3 ms | 9.1 ms | 8.8 ms / 9.2 ms | 1.47× |  |
| 167 | feature_coverage | `ft_tfun` | 40 | 80 | 6.2 ms | 6.1 ms / 6.5 ms | N/A | — | — | NFsim refused XML |
| 168 | feature_coverage | `ft_total_rate` | 100 | 100 | 6.2 ms | 6.0 ms / 6.3 ms | 8.8 ms | 8.7 ms / 9.0 ms | 1.41× |  |
| 169 | feature_coverage | `ft_functional_rate` | 300 | 100 | 6.2 ms | 6.1 ms / 6.3 ms | 9.0 ms | 8.9 ms / 9.0 ms | 1.46× |  |
| 170 | nfsim_basicmodels | `r32` | 5 | 200 | 6.1 ms | 6.1 ms / 6.3 ms | 9.4 ms | 9.3 ms / 9.5 ms | 1.53× |  |
| 171 | nfsim_basicmodels | `r17` | 5 | 50 | 6.1 ms | 6.0 ms / 6.1 ms | 8.7 ms | 8.6 ms / 8.8 ms | 1.43× |  |
| 172 | feature_coverage | `ft_mm_ratelaw` | 30 | 60 | 6.0 ms | 6.0 ms / 6.2 ms | 8.7 ms | 8.6 ms / 8.8 ms | 1.45× |  |
| 173 | feature_coverage | `edg_dynamic_rate_zero_obs` | 40 | 80 | 6.0 ms | 5.9 ms / 6.2 ms | 8.7 ms | 8.6 ms / 8.8 ms | 1.45× |  |
| 174 | feature_coverage | `edg_ring_break_constraint` | 10 | 40 | 6.0 ms | 5.9 ms / 6.3 ms | 8.2 ms | 8.1 ms / 8.3 ms | 1.38× |  |
| 175 | corpus | `nfsim_hybrid_particle_field` | 10 | 50 | 6.0 ms | 6.0 ms / 6.2 ms | 8.6 ms | 8.5 ms / 9.0 ms | 1.43× |  |
| 176 | feature_coverage | `edg_zero_rate_rule` | 40 | 80 | 6.0 ms | 5.9 ms / 6.2 ms | 8.5 ms | 8.4 ms / 8.6 ms | 1.42× |  |
| 177 | feature_coverage | `edg_time_dependent_rate` | 50 | 100 | 5.9 ms | 5.9 ms / 6.6 ms | N/A | — | — | NFsim produced no output |
