# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

- **Multi-mol Molecules observable counts on palindromic patterns
  with symmetric components** (`5d90724`). The BFS in
  `count_multi_molecule_embeddings` committed to the first valid
  partner embedding consistent with the walked bond and silently
  dropped sym-equivalent alternatives. On a 5-mol palindromic
  observable like `B(c!1).R(b!1,a!2).A(r!2,r!3).R(b!4,a!3).B(c!4)`
  with a sym partner reached via a non-sym bond endpoint (the
  basicmodels v07 / r07 shape), this under-counted by exactly the
  partner's symmetry factor. Replaced the BFS body with a recursive
  enumerator that branches over every partner embedding consistent
  with the walked bond. Untouched: rule rates (different code path),
  `Species` observables (deduped by complex anyway), single-mol
  observables, and the 2-mol-1-bond fast path.

- **Multi-mol unimolecular rule rates over-counted on hosts with
  symmetric components** (`7472a07`). `count_multi_mol_fast` (the
  generic path; the 2-mol-1-bond specialization is gated off for
  this shape) called `count_embeddings_single` for the seed without
  `reacting_local`, so injective embeddings differing only in
  non-reacting sym slots were never deduped. On the basicmodels
  v02 / r02 unbind rule `X(y!1, p~0).Y(x!1) -> X(y, p~0) + Y(x)`
  this fired at 2× the BNG2-strict rate. Threaded `reacting_local`
  through `count_multi_mol_fast` → `count_multi_mol_fast_generic` →
  the seed-side `count_embeddings_single` call. Phos remains
  correct (mult=2): StateChange targets the matched p, so the two
  sym embeddings produce different keys under dedup and stay
  distinct.

- **Compile-time embedding correction was double-applied for
  multi-mol patterns after the seed-dedup fix** (`3423b0d`). The
  previous fix subsumed the work `compute_embedding_correction_multimol`
  was doing, so applying both halved the rate when the pattern had
  sym non-reacting components (the basicmodels v18 / r18 shape).
  Set `embedding_correction_a/_b = 1.0` for multi-mol (mirroring
  single-mol, which always did this). Removed
  `compute_embedding_correction[_multimol]` and `_impl` —
  ~110 lines of dead code.

### Removed

- **Four upstream NFsim regression tests dropped from the basicmodels
  suite** (`85feae1`, `9fb2efb`). All four tested NFsim-specific
  behaviors that don't apply to RuleMonkey:
  - `r33` and `r35` pinned NFsim issues #22 / #21 ("occupied-site
    bond error") and #14 (RHS `.` between products that NFsim
    splits anyway). On the BNGL these tests carry, BNG2.pl's
    `generate_network` produces the chemistry-correct behavior
    (bound by free-B count for r33; zero reactions for r35) and
    RuleMonkey matches BNG2.pl. The NFsim references captured the
    historic NFsim quirks, which by design diverge from BNGL strict.
  - `r31` is a crash regression test with no `begin observables`
    block (the author's own comment: *"validation harness will run
    NFsim on this XML and ensure it doesn't crash"*).
  - `r34` includes a `begin observables` block but the author
    deliberately commented out the only line in it.

  After the removals the suite is clean **29 PASS / 0 FAIL /
  0 NO_MATCH**. Joins the pre-existing r27 / r28 / r36 in the
  PROVENANCE appendix.

### Added

- Three feature_coverage regression tests pinning the sym-K shapes
  fixed above:
  - `ft_multimol_sym_obs.bngl` — 5-mol palindromic observable with
    sym partner (the r07 shape).
  - `ft_multimol_unimol_unbind_sym.bngl` — multi-mol unimolecular
    unbind on a host with sym components, where operations don't
    differentiate the embeddings (the r02 shape).
  - `ft_multimol_pattern_sym_nonreacting.bngl` — multi-mol pattern
    with sym non-reacting components on the seed (the r18 shape,
    catches the embedding-correction × dedup double-apply).
  Each was verified to fail pre-fix and pass post-fix via stash-
  and-rerun.

### Changed

- `tests/reference/basicmodels/PROVENANCE.md` rewritten to lead with
  what the suite *is* (29 imported tests, source, reference flow)
  and treat the seven upstream NFsim tests not carried over as a
  clearly-labeled appendix grouped by reason. `harness/basicmodels.py`
  and `harness/generate_basicmodels_refs.py` docstrings tightened
  to a one-line origin pointer.

## [3.0.0] — 2026-04-26

First release of the cleanroom C++17 rewrite. Not source-, ABI-, or
CLI-compatible with RuleMonkey 2.0.25.

### Added
- Cleanroom C++17 simulation engine (`librulemonkey`) with no external
  dependencies beyond the standard library — including an in-house XML
  parser, so neither TinyXML nor any other third-party XML library is
  required.
- Public C++ API at `include/rulemonkey/{simulator,types}.hpp` exposing
  `rulemonkey::RuleMonkeySimulator` for in-process embedding (init / run /
  step_to / simulate / add_molecules / set_param / save_state / load_state).
- `rm_driver` batch CLI emitting `.gdat`-format trajectories.
- Test corpora under `tests/models/`:
  - `feature_coverage/` — 51 BNGL feature-coverage models with invariants
    and golden values.
  - `corpus/` — 71 real-world rule-based models for efficiency and
    correctness benchmarking.
  - `nfsim_basicmodels/` — 72 models from the NFsim parity suite.
- 100-replicate ensemble reference trajectories under
  `tests/reference/nfsim/` (regenerated 2026-04-10) — 69 of 71 models
  from a gold-standard NFsim build, with `toy_jim` and `rm_tlbr_rings`
  replaced by hand-rolled Gillespie SSA where NFsim was confirmed to
  produce incorrect output. SSA scripts under `harness/ssa/`. See
  `tests/reference/nfsim/PROVENANCE.md` for the full provenance and
  per-model exceptions table.
- Python harness scripts under `harness/` for end-to-end benchmarking
  and validation.
- CMake (≥3.25) build with Ninja generator and configurable presets;
  smoke test wired into CTest.
- GitHub Actions CI for Linux and macOS.

### Changed
- License: **GPLv3 → MIT.** Cleanroom code reuses no legacy source.
- Build system: **autotools → CMake** with Ninja.
- Engine language: **C → C++17.**
- Default semantics: strict BNGL `block_same_complex_binding` is now on
  by default; pass `-no-bscb` for compatibility with NFsim runs that
  omitted `-bscb`.

### Removed
- Legacy C implementation (RuleMonkey 2.0.25). Preserved in the parent
  fork's git history (`richardposner/RuleMonkey`).
- Vendored `dSFMT-1.3` and `nauty22`. The cleanroom uses
  `std::mersenne_twister_engine` for RNG and **does not currently
  exploit graph canonical labeling** for complex/species identification.
  Re-introducing canonical labeling is a candidate optimization for a
  future minor release.

### Known limitations
- No Python bindings yet. Python access is via the `rm_driver`
  subprocess. Native bindings (`pybind11` + `scikit-build-core`) are
  planned for 3.1.0 alongside the BNGsim integration.
- Compartments are refused at Tier 0 (exit code 2). Volume scaling and
  surface chemistry remain open work.
- Arrhenius rate laws and a small set of unsupported BNGL constructs
  cause hard refusals; pass `--ignore-unsupported` to demote to
  warnings for testing.

### Lineage
The legacy implementation, RuleMonkey 2.0.25, was introduced in:

> Colvin J, Monine MI, Gutenkunst RN, Hlavacek WS, Von Hoff DD, Posner
> RG. *RuleMonkey: software for stochastic simulation of rule-based
> models.* BMC Bioinformatics 11:404 (2010). PMID: 20673321.

[3.0.0]: https://github.com/wshlavacek/RuleMonkey/releases/tag/v3.0.0
