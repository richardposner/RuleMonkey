# basicmodels Reference Data — Provenance

100-replicate NFsim ensemble references for the historical NFsim parity
test suite (33 testable models). Three models from the original NFsim
test set were removed because they don't apply to RM:

- **r27, r28** — used the BNGL `population` keyword (hybrid
  particle-population SSA, Hogg 2013). RM has no equivalent and now
  refuses such models at Tier-0 with a `MoleculeType@population`
  error pointing at the relevant `<MoleculeType>` element. See
  `cpp/rulemonkey/simulator.cpp:scan_unsupported`.
- **r36** — tested NFsim's `-gml auto` fallback (issue #17). RM only
  honors numeric `set_molecule_limit`; the `auto` token is a
  non-feature, not a bug.

The original `r{27,28,36}.txt` and `v{27,28,36}.bngl` source files
were deleted from `tests/models/nfsim_basicmodels/`.

**Last regeneration:** 2026-04-26
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
