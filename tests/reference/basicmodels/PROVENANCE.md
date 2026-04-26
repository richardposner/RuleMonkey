# basicmodels Reference Data — Provenance

100-replicate NFsim ensemble references for the historical NFsim parity
suite (33 of 36 models; 3 excluded — see below).

**Last regeneration:** 2026-04-26
**NFsim binary:** `/Users/wish/Code/nfsim-rm/build/NFsim` (has UTL+1 fix)
**BNG2:** `/Users/wish/Simulations/BioNetGen-2.9.3/BNG2.pl`
**Reps per model:** 100

## Excluded models

- **r27** — uses BNGL `population` keyword (hybrid particle-population SSA; not in RM)
- **r28** — uses BNGL `population` keyword (hybrid particle-population SSA; not in RM)
- **r36** — tests `-gml auto` (NFsim issue #17 fallback for large initial pops; not in RM)

## Flag translation policy

The `r{NN}.txt` config files specify NFsim run options that include a
mix of trajectory-affecting flags and output-only / debug flags. The
generator translates them so the gold-standard reference reflects what
the trajectory should look like, not which output files NFsim happens
to produce. See `harness/generate_basicmodels_refs.py:STRIP_FLAGS` for
the canonical list.

Notably, `-utl 3` (used by r21–r26 in the original config) is dropped
on both sides — both NFsim (with the UTL+1 patch) and RM auto-compute
`max_pattern_size + 1`, which is the correct value.

## Per-model summary

See `sim_params.tsv` for the actual flags used per model and the wall
time of each NFsim ref-gen run.
