# basicmodels Reference Data — Provenance

NFsim-generated 100-replicate ensemble references for the
`tests/models/nfsim_basicmodels/` BNGL test suite. RuleMonkey runs each
model and verifies its trajectories are statistically indistinguishable
from the NFsim ensemble (per-rep z-score over the trapezoidal integral
of each observable).

## Source

The BNGL files in `tests/models/nfsim_basicmodels/` are derived from
NFsim's own `validate/basicModels/` regression suite — small,
purpose-built tests that NFsim used to verify specific behaviors of its
rule engine. The suite was inherited when RuleMonkey was forked from
NFsim and is preserved here as one of three test corpora (alongside
`feature_coverage/` and the larger model corpus under `corpus/`).

The 29 BNGL files cover: simple binding / unbinding / phosphorylation,
multiple identical sites and dimerization, ring closure (single- and
two-molecule-type), explicit and implicit cross-phosphorylation,
species deletion / synthesis, mRNA monomer connectivity (r21–r26),
function and parameter handling, and a few intramolecular and
edge-case feature checks (r29, r30, r32). Per-model names, simulation
times and flags are in `sim_params.tsv`.

## Reference generation

`harness/generate_basicmodels_refs.py` per model:
  1. Runs BNG2.pl on `v{NN}.bngl` and takes the resulting XML.
  2. Translates the NFsim flags from `r{NN}.txt` — strips output-only
     flags and the explicit `-utl 3` override so both engines use the
     auto-computed `max_pattern_size + 1` (the correct value).
  3. Runs NFsim 100 times with seeds 1..100.
  4. Aggregates to `mean.tsv` and `std.tsv` per observable per time
     point. Per-rep `.gdat` files land in `replicates/r{NN}/`
     (gitignored).

**NFsim binary:** a local NFsim build patched with the UTL+1 fix
(NFsim's auto-computed unbinding-trace-length is one less than the
correct `max_pattern_size + 1`; without the patch, large-pattern
basicmodels rules silently miss matches).  Set the `NFSIM_BIN` env
var when regenerating; an unpatched binary will produce diverging
references for r07-style multi-mol patterns.
**BNG2:** `~/Simulations/BioNetGen-2.9.3/BNG2.pl`.
**Reps per model:** 100.
**Last regenerated:** 2026-04-26 (full 33-model run); r31, r33, r34,
r35 were later dropped on 2026-04-27 — see appendix.

## Flag translation policy

NFsim's `r{NN}.txt` configs mix trajectory-affecting flags with
output-only / debug flags. The generator translates them so the
gold-standard reference reflects what the trajectory should look like,
not which output files NFsim happens to produce. See
`harness/generate_basicmodels_refs.py:STRIP_FLAGS` for the canonical
list.

## Per-model summary

See `sim_params.tsv` for the per-model flags actually used and the
wall time of each NFsim ref-gen run.

---

## Appendix — NFsim upstream tests not included

NFsim's `validate/basicModels/` originally has 36 BNGL files. Seven
aren't carried over here because they test features or behaviors that
don't apply to RuleMonkey. They aren't failures and they aren't pending
work — there is simply nothing useful to verify on the RM side, so
they were never imported / were removed once the situation was clear.

Grouped by reason:

**Features RM doesn't implement**
- `r27`, `r28` — the BNGL `population` keyword (hybrid
  particle-population SSA from Hogg 2013). RuleMonkey is a pure
  network-free SSA and refuses such models at Tier-0 with a
  `MoleculeType@population` error
  (`cpp/rulemonkey/simulator.cpp:scan_unsupported`).

**NFsim flags RM doesn't expose**
- `r36` — exercises NFsim's `-gml auto` flag (NFsim issue #17
  fallback). RuleMonkey only honors the numeric
  `set_molecule_limit(N)` form, so there is no equivalent path to
  exercise.

**NFsim regression tests for bugs that contradict BNGL strict semantics**
- `r33` — pinned NFsim issue #22 / #21 ("occupied-site bond error").
  The rule
  `B(Auto) + C(sub!1).A(ub!1) -> B(Auto!1).C(sub!1) + A(ub)` requires
  a free `B(Auto)`, but NFsim continues firing after the only B
  becomes bonded — applying the DeleteBond half (breaking A-C) while
  silently dropping the AddBond half. BNG2.pl `generate_network`
  correctly bounds the rule by the free-B pool, and RuleMonkey
  matches BNG2.
- `r35` — pinned NFsim issue #14. RHS uses `.` between products to
  indicate "same complex" (e.g. `B(b1!1).C(c1!1) -> B(b1).C(c1)`),
  but NFsim ignores the constraint and splits anyway. BNG2.pl
  `generate_network` emits zero reactions for that rule shape
  because the products would not be in one complex; RuleMonkey
  matches that strict behavior (the rule fires only when the bond
  is in a cycle and the dot semantics can be honored).

The NFsim references for `r33` and `r35` therefore record the
historical NFsim behavior, which by design diverges from BNGL strict.
Including them in the parity suite would test that RuleMonkey
reproduces NFsim bugs, which is not a goal.

**Tests with no comparable observables**
- `r31` — a crash-regression test. The author's own comment is
  explicit: *"validation harness will run NFsim on this XML and ensure
  it doesn't crash."* The model has no `begin observables` block and
  was never intended to produce a comparable trajectory. Crash /
  robustness coverage is the job of `ctest` and the smoke harness,
  not the parity suite.
- `r34` — pinned NFsim issue #24 ("RandomBiBi occupied-site
  edge case"). The author included a `begin observables` block but
  deliberately commented out the only line in it (`# Molecules Q
  Q()`), so there is no observable to integrate against. Like r33
  and r35, it's a regression test for an NFsim-specific edge case
  rather than a chemistry test that would generalize across
  engines.
