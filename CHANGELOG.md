# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
  `tests/reference/nfsim/` from a gold-standard NFsim build (regenerated
  2026-04-10; see `tests/reference/nfsim/PROVENANCE.md`).
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
