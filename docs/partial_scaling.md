# Partial scaling (approximate acceleration)

Partial scaling is an **opt-in, approximate** simulation mode that trades a
controlled amount of statistical accuracy for fewer SSA steps. It is **off by
default** — an unset critical population runs the exact network-free SSA, bit
-for-bit identical to a build without the feature.

> **Approximate.** With partial scaling on, individual trajectories are *not*
> samples from the exact stochastic process. Ensemble **first moments (means)
> stay unbiased**, but **variance is inflated** by an amount that shrinks as the
> critical population grows (`≈ C/Nc`). Do not use it where you need exact
> single-trajectory statistics or exact higher moments. See
> [Accuracy](#accuracy-what-you-trade-away) below.

## Enabling it

CLI (`rm_driver`):

```bash
# --nc <N> : critical population.  N<=0 (default) = exact SSA.
build/release/rm_driver model.xml 10 100 7 --nc 50 > traj.gdat
```

In-process API:

```cpp
rulemonkey::RuleMonkeySimulator sim("model.xml");
sim.set_critical_population(50);   // approximate mode; throws if a session is active
auto r = sim.run({0.0, 10.0, 100}, /*seed=*/7);
```

`set_critical_population` mirrors `set_molecule_limit`: instance-local, applied at
the next `run()` / `initialize()`, and it throws if a session is already active.
It is a no-op (`Nc <= 0`) unless you set a positive value.

## What it does

An exact network-free SSA fires **one** reaction per step. With a critical
population `Nc`, a rule whose reactant population `n_r` exceeds `Nc` fires a
**batch** of `K_r = ⌈n_r / Nc⌉` reactions under a **single** `incremental_update`
of the affected state (particle batching / "Level-A" partial scaling). The
rule's propensity is scaled by `λ_r = 1/K_r` so the total propensity — and the
time step `dt` — stay consistent with the exact process's first moment.

Concretely, a rule at population `n_r = 5000` with `Nc = 50` fires in batches of
`K_r = 100` reactions per step. A rule whose population is at or below `Nc`
keeps `K_r = 1` and fires exactly, one reaction at a time. So `Nc` is the
population **below which a rule is treated exactly** — smaller `Nc` ⇒ larger
batches ⇒ more aggressive acceleration ⇒ more excess variance.

The per-rule, time-averaged batch multiplier `m_r = Σ K_r·Δt / Σ Δt` is reported
on `stderr` when `--nc` is set (and on `Result::ps_multipliers`), so you can see
which rules actually batched:

```
[ps] nc=50 firings=27064 steps=11771 mult=2.5,11.7,11.7,21.1,...
```

`firings` is the number of reactions applied; `steps` is the number of SSA
iterations. `firings/steps` is the mean batch size; `firings > steps` is the
positive signal that batching happened.

## What it buys

Partial scaling reliably reduces the **SSA step count** — each batched step does
the work of `K` reactions under one update. It does **not** reduce the *total*
reaction work: the batched `incremental_update` runs over the union of all `K`
firings' affected molecules, so batching **repackages** the per-reaction update
cost rather than removing it. The wall-clock win is therefore the amortized
per-*event* fixed cost (rule selection, RNG, epoch bookkeeping) spread over `K`
firings — real, but modest on RuleMonkey's already-optimized C++ engine, and
model-dependent.

Measured on the RuleMonkey benchmark models at an aggressive critical
population (`Nc = 5`; `Nc = 10` for the two ~500-copy ring toys), 5 reps each,
`rm_driver --release`. **SSA-step reduction** = exact steps / scaled steps;
**wall speedup** = exact wall / scaled wall (min of reps, to reject transient
machine load):

| model         | SSA-step reduction | wall speedup |
|---------------|-------------------:|-------------:|
| `lat`         |               347× |        1.21× |
| `homodimer`   |                50× |        0.97× |
| `intra_ring`  |                49× |        1.12× |
| `ab_ring`     |                40× |        0.95× |
| `machine`     |                34× |        1.00× |
| `tcr`         |                28× |        1.21× |
| `egfr_net`    |               9.9× |        1.01× |
| `birth_death` |               1.9× |        1.04× |

The pattern: **large step-count reduction, roughly wall-neutral (0.85–1.21×).**
Models with many cheap rules (e.g. `tcr`) benefit most from the fixed-cost
amortization; complex-heavy models with large per-reaction updates (e.g. `lat`,
`machine`) see the update work dominate, so their wall time is nearly flat — and
can even dip *below* 1× at moderate `Nc`, where the bigger union update outweighs
the fixed-cost saving. The diagnostic is how much a batched step costs versus the
firings it replaces: on `tcr` at `Nc=5` a batched step costs ~23× an exact step
while doing the work of ~28 firings, and that gap is the fixed-cost saving you
feel as the 1.21× speedup; on `machine` a batched step costs ~90× for ~96
firings — essentially all repackaged, hence ~1.00×.

> The Python `nfa` prototype this feature is ported from shows a *large*
> wall-clock speedup, because an interpreted engine's per-event fixed cost is
> enormous. RuleMonkey has already minimized that cost, so on RuleMonkey the
> headline payoff is step-count reduction (and the foundation it lays for
> future population-scaling work), not raw wall time.

## Accuracy: what you trade away

- **First moment (mean) — unbiased.** Ensemble means of every observable match
  the exact-SSA ensemble (and the analytic equilibrium where one exists) within
  statistical error at every tested `Nc`, down to `Nc = 5` on large-population
  models. (On a *small*-population observable — e.g. a species with ~10 copies —
  aggressive `Nc` introduces a mild finite-population bias, inherent to the
  method and shared with `nfa`.)
- **Second moment (variance) — inflated, shrinking in `Nc`.** Partial scaling
  adds excess variance `≈ C/Nc`. It is largest at small `Nc` and decreases
  monotonically toward the exact variance as `Nc` grows. For `birth_death`
  (`A ~ Poisson(100)`, exact variance ≈ 100; means unbiased throughout,
  `|z| < 1.5`), across 300 reps:

  | `Nc` | variance | inflation vs exact |
  |------|---------:|-------------------:|
  | exact (off) |  106 | 1.0× |
  | 50   |      131 | 1.2× |
  | 20   |      290 | 2.7× |
  | 10   |      550 | 5.2× |
  | 5    |     1116 | 10.6× |

  The inflation scales with the **batched population**: a given `Nc` is gentler
  on a 100-copy species than on a 1000-copy one (a 500-copy homodimer at `Nc=50`
  inflates variance ~11×, the same factor `birth_death` reaches only at `Nc=5`).
  Because `Nc` is an absolute population, always read it relative to your model's
  species counts.

## Choosing `Nc`

There is no universal best value — it depends on your population sizes and how
much variance inflation you can tolerate. Guidance:

- **`Nc` is the exact-below-this population.** Set it well below your smallest
  *dynamically important* species population so those species stay exact, but
  small enough that your large populations batch.
- **Speed saturates; variance grows.** Wall-clock speedup rises as `Nc` shrinks
  but saturates quickly (bounded by the fixed-cost fraction), while excess
  variance keeps growing as `1/Nc`. There is little reason to push `Nc` below
  the knee.
- **Default recommendation: `Nc ≈ 20–50`** for populations in the hundreds. It
  batches large populations several-fold — capturing most of the (modest) wall
  benefit — while holding variance inflation to roughly 1.5–3× and staying clear
  of the guard on all but the slowest models. Drop to `Nc ≈ 5–10` only when you
  want maximum step-count reduction and can tolerate ~5–10× variance inflation;
  raise it when accuracy matters more than speed. Scale the number to your
  model's populations (see the variance note above) — there is no universal
  constant.
- Validate on your own model: run a small ensemble at your chosen `Nc` and at
  `Nc` off (exact) and confirm the means agree and the variance inflation is
  acceptable before trusting a production sweep.

## The `Nc`-too-small guard

When scaling is active, the total propensity is scaled down by `λ_r`. If `Nc` is
so small (batches so large) that the expected time to the next reaction already
spans the whole remaining window, the run would collapse to a single step and
produce garbage. Rather than silently return nonsense, the run **fails loudly**:

```
ERROR: Partial scaling: Nc=5 is too small for this model — the scaled total
propensity 0.0189 implies an expected step dt=52.9s, which meets or exceeds the
remaining window 10s (t_end=10). The run would collapse to a single step.
Increase Nc.
```

This is a one-time check at the start of the run. It fires on **slow-dynamics /
low-propensity models with a short time window** (e.g. gelling ring models like
`rm_tlbr_rings`), where aggressive scaling is inappropriate — increase `Nc` (or
lengthen `t_end`). Fast models with high propensity do not trip it even at small
`Nc`. The guard only runs when scaling is actually active; a huge `Nc` (an exact
run) and a model whose own dynamics are naturally slow never trip it.

## Not an oracle for parity

NFsim does **not** implement this partial-scaling method, so an NFsim run is not
a reference for a scaled RuleMonkey run. The scaled path is validated against
exact RuleMonkey, analytic moments, the `nfa` prototype, and BNG2 SSA/ODE
ensembles — never against NFsim. With `Nc` unset, RuleMonkey remains bit-for-bit
identical to NFsim-parity exact behavior.
