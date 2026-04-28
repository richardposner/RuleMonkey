# basicmodels Reference Data — Provenance

100-replicate NFsim ensemble references for the historical NFsim parity
test suite (31 testable models). Five models from the original NFsim
test set were removed — three for missing-feature reasons, and two
because their reference data captures NFsim-specific quirks that
contradict BNG2.pl's strict BNGL semantics (RM matches BNG2.pl):

- **r27, r28** — used the BNGL `population` keyword (hybrid
  particle-population SSA, Hogg 2013). RM has no equivalent and now
  refuses such models at Tier-0 with a `MoleculeType@population`
  error pointing at the relevant `<MoleculeType>` element. See
  `cpp/rulemonkey/simulator.cpp:scan_unsupported`.
- **r36** — tested NFsim's `-gml auto` fallback (issue #17). RM only
  honors numeric `set_molecule_limit`; the `auto` token is a
  non-feature, not a bug.
- **r33** — pinned NFsim issue #22 / #21 ("occupied-site bond
  error"). The rule
  `B(Auto) + C(sub!1).A(ub!1) -> B(Auto!1).C(sub!1) + A(ub)` requires
  a free `B(Auto)`, but NFsim continues firing after the only B
  becomes bonded — running the DeleteBond half (breaking A-C) while
  silently dropping the AddBond half. BNG2.pl's `generate_network`
  correctly bounds the rule by the free-B pool, and RM matches that
  bound. The NFsim reference captures the buggy unbounded firing.
- **r35** — pinned NFsim issue #14 (rule RHS uses `.` between
  products to indicate "same complex" but NFsim ignores it and
  splits). On `B(b1!1).C(c1!1) -> B(b1).C(c1) kr`, BNG2.pl
  `generate_network` emits zero reactions because the product `.`
  pattern requires the two molecules to remain in one complex,
  which is impossible after deleting their only bond. RM matches
  BNG2.pl's strict refusal (the rule fires only when the bond is
  in a cycle, e.g. on the original triangle); NFsim fires
  unconditionally.

The original `r{27,28,33,35,36}.txt` and `v{27,28,33,35,36}.bngl`
source files were deleted from `tests/models/nfsim_basicmodels/`.

**Last regeneration:** 2026-04-27 (r33/r35 entries removed; rest
unchanged from 2026-04-26)
**NFsim binary:** `/Users/wish/Code/nfsim-rm/build/NFsim` (has UTL+1 fix)
**BNG2:** `/Users/wish/Simulations/BioNetGen-2.9.3/BNG2.pl`
**Reps per model:** 100

## Flag translation policy

The `r{NN}.txt` config files specify NFsim run options that include a
mix of trajectory-affecting flags and output-only / debug flags. The
generator translates them so the gold-standard reference reflects what
the trajectory should look like, not which output files NFsim happens
to produce. See `harness/generate_basicmodels_refs.py:STRIP_FLAGS` for
the canonical list.

Notably, `-utl 3` (used by r21–r26 in the original config) is dropped
on both sides — both NFsim (with the UTL+1 patch from
`nfsim-rm` commit `40e6b93`) and RM auto-compute
`max_pattern_size + 1`, which is the correct value.

## Per-model summary

See `sim_params.tsv` for the actual flags used per model and the wall
time of each NFsim ref-gen run.
