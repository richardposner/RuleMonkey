# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Removed

- **Dead pytest config in `pyproject.toml`.** The
  `[tool.pytest.ini_options]` block declared
  `testpaths = ["tests", "harness"]`, but `tests/` is C++-only and
  `harness/` holds benchmark/research scripts — no `test_*.py`,
  `*_test.py`, or `conftest.py` exists anywhere.  Running `pytest`
  from the repo root collected zero tests, which is a misleading
  signal for an external reader.  Removed the config block and the
  unused `pytest>=8.0` dev dependency; the suite is C++ ctest plus
  harness-driven Python scripts.

- **Dead helper functions and locals in `simulator.cpp`.**
  Three unused free functions in the XML-parser anonymous namespace
  (`need_child`, `has_attr`, `any_rule_has_child`) and two unused
  local variables (`rp_start_0`, `rp_start_1` inside the
  same-components detection) had been triggering compiler
  `-Wunused-function` / `-Wunused-variable` warnings on every clean
  build.  Removed; the library now compiles warning-free under both
  the `release` and `asan` presets.

### Changed

- **Parameter forward-reference resolution iterates to fixed point.**
  `load_model` previously did a single retry pass after the initial
  resolution pass, which handled at most one level of forward
  reference.  BNG2 emits parameters in dependency order so this
  never bit in practice, but arbitrary XML declaring
  `P3 = 2*P2; P2 = P1; P1 = 1` in that order would have left `P3`
  stale.  Now iterates until either every value is stable or the cap
  (param-count + 4) is hit.

- **CONTRACT comment on `*_expr` fields in `model.hpp`.**  Each
  symbolic-source field (`rate_expr`, `mm_kcat_expr`, `mm_Km_expr`,
  `concentration_expr`) now carries an inline note saying it must be
  re-resolved in `RuleMonkeySimulator::Impl::apply_overrides`.
  Defensive against future parser extensions that add a 5th
  parameter-derived field and forget to wire it into the override
  cascade — exactly the regression that the 3.1.x `apply_overrides`
  fix was designed to close.

- **`asan` preset disables `RULEMONKEY_INSTALL`; CMakeLists refuses
  `RULEMONKEY_INSTALL=ON` with `RULEMONKEY_ENABLE_ASAN=ON`.**  The asan
  flags are PRIVATE compile / PUBLIC link on the `rulemonkey` target,
  so an asan-instrumented installed archive would propagate
  `-fsanitize=address,undefined` to every downstream
  `find_package(RuleMonkey)` consumer's final link line — failing
  outright unless the consumer also enables asan, and producing a
  binary that runs against the wrong runtime even when it links.
  Asan is dev-only; the supported `asan` preset now auto-disables
  install, and the build refuses the dangerous combo with an
  actionable error message pointing at `-DRULEMONKEY_INSTALL=OFF`.

### Fixed (docs)

- **CHANGELOG 3.0.0 model count off by one.** The 3.0.0 entry said
  "51 BNGL feature-coverage models"; `git ls-tree v3.0.0 --
  tests/models/feature_coverage/` counts 52 `.bngl` files.  Corrected.

- **README pointer to the find_package / add_subdirectory snippet.**
  The "Embedding (C++ API)" section pointed readers at
  `examples/CMakeLists.txt` for the CMake consumption snippet, but
  that file is two lines (`add_executable` + `target_link_libraries`).
  The actual snippets live in the doc-comment header of
  `examples/embed.cpp`.  Updated the pointer.

- **`docs/model_semantics.md` "RM v1" wording.**  The Tier-0 refusals
  table for multi-molecule and duplicate Fixed species said "RM v1
  only supports …" / "RM v1 allows at most …".  RM is at 3.x; a
  reader could mistake "v1" for RuleMonkey 1.x rather than the v1
  scope of the Fixed-species feature itself.  Wording tightened to
  "RM currently …".

### Fixed

- **`get_parameter` returned the parsed-at-load value between
  `set_param` and the next `run()` / `initialize()`.** The
  override map only synced into `model.parameters` inside
  `apply_overrides()`, which only ran on session start. Embedders
  querying overrides via `get_parameter()` between runs saw stale
  values. `set_param` and `clear_param_overrides` now invoke a
  light `sync_parameters()` so the public read is coherent
  immediately.

- **Derived parameter expressions (`B = 2*A`) did not cascade
  through `set_param` overrides.** Parameter values were resolved
  once at XML parse time and the override map only splatted the
  literal name's value, leaving downstream derivations frozen at
  their parsed numeric. Captured the symbolic expression for each
  declared parameter at parse time (`Model::parameter_exprs`) and
  re-cascade in declaration order inside `sync_parameters()` so
  set_param on a base parameter propagates to every derived
  parameter that references it. Override on a derived parameter
  still wins (skips the expression for that name).

### Changed

- **`set_param` validates parameter names against the loaded XML.**
  Typo'd names previously leaked into `param_overrides` as silent
  no-ops; they now throw `std::runtime_error("Unknown parameter
  '...'")`. Mirrors `get_parameter`'s existing throw on unknown
  names.

- **`docs/model_semantics.md`**: new "Parameter overrides" section
  documents the cascade behavior and the unknown-name rejection.

### Added

- Regression coverage in `tests/cpp/set_param_test.cpp` and a new
  `tests/cpp/derived_param_model.xml` fixture — two-layer derived-
  parameter chains (`A_tot = A_base * A_factor`,
  `kp = kp_base * kp_mult`) plus get_parameter-coherence and
  unknown-name-rejection cases. Each new assertion was verified to
  fail on the unfixed engine via stash-and-rerun.

## [3.1.0] — 2026-04-29

### Fixed

- **Self-binding propensity under bscb under-counted by `(N-1)/N`**
  (`9f25bba`). `compute_propensity` was subtracting within-molecule
  pair contributions (`extra_eff`) at propensity time AND
  `select_reactants` was rejecting `mol_a == mol_b` at sample time —
  the same self-pairs were removed twice. Hetero-binding was
  unaffected (`extra_eff` is naturally zero across distinct molecule
  types); the bias only appeared on `A+A` shapes (both symmetric
  `A(c)+A(c)` and asymmetric same-type `A(l)+A(r)`). Hidden under
  the standard 20-rep 5σ benchmark noise floor; emerged at 100+
  reps. RM-vs-NFsim avg-z dropped from ~−1.0 to ~−0.2 across the
  affected models post-fix; vs deterministic ODE on
  `combo_addbond_connected`, RM moved from z=−3.21 to z=−1.64 (now
  better than NFsim's z=−2.02). Surfaced by the new `edg_*` stress
  suite at 500 reps. Same bias also affected the existing
  `combo_addbond_connected` corpus model.

- **MM(kcat, Km) rate law was silently zero-rate** (`6421017`).
  The XML loader parsed `RateLawType::MM` and stored `mm_kcat` /
  `mm_Km` on the rule, but the engine never read those fields —
  every MM rule ran at zero propensity. Confirmed live: an MM rule
  on a 100-substrate / 5-enzyme system left S unchanged at 100
  forever in RM, while NFsim depleted S to ~10 by t=30. Now uses
  NFsim's QSS formula `sFree = 0.5·((S−Km−E) + √((S−Km−E)² + 4·Km·S))`,
  `a = kcat·sFree·E/(Km+sFree)` (mirrors `MMRxnClass::update_a` in
  the NFsim source). Verified at 50 reps each: RM and NFsim agree
  within stochastic noise.

- **TFUN sentinel substitution + .tfun file resolution**
  (`cb03071`). Two bugs in RM's TFUN code path, both surfaced by
  the new `ft_tfun.bngl` test:
  - Sentinel name mismatch: RM looked for `__TFUN_VAL__` (single
    underscore) but BNG2 2.9.3 emits `__TFUN__VAL__` (double
    underscore). The dev branch `fix-tfun-has-tfuns-reset` reverts
    to the single-underscore form for the new lowercase `tfun()`
    syntax. RM now accepts both, longer match first.
  - File-resolution path: RM only searched relative to the XML
    directory; the harness lays out XML in `tests/.../xml/` with
    the `.tfun` next to the source BNGL one level up. RM now falls
    back to `<xml_dir>/..` so author-side and harness-side
    layouts both work.

- **Multi-mol Molecules observable counts on palindromic patterns
  with symmetric components** (`5d90724`). The BFS in
  `count_multi_molecule_embeddings` committed to the first valid
  partner embedding consistent with the walked bond and silently
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

### Added

- **`edg_*` stress-test suite — 20 honest probes designed to break
  RM** (`38d2087`, `6617a84`, `563495b`). Targets feature
  combinations not covered by the existing `ft_*` / `combo_*`
  suite: state-increment ladders, synthesis-into-pre-bonded
  complexes, time-explicit rates, pattern-level local functions,
  ternary embeddings, branched aggregates, multi-Fixed competition,
  zero-rate edge cases, and several `A+A` self-binding shapes that
  were hiding the propensity bug above. 18 models use BNG2 ODE as
  the deterministic verdict reference; the two polymer-style models
  (`edg_pattern_local_fcn`, `edg_branched_polymer`) use 100-rep
  NFsim refs since their network generation is intractable. Per-
  model rationale and 500-rep verification table at
  `tests/models/feature_coverage/EDG_RATIONALE.md`.

- **Tier-0 refusal for Sat / Hill / FunctionProduct rate laws**
  (`13bb424`). The rule loader recognised only `Ele`, `Function`,
  and `MM`; everything else fell through to the default Ele type
  with `rate_value=0.0`, so rules using `Sat()` / `Hill()` /
  `FunctionProduct()` silently produced zero propensity and never
  fired. RM now refuses loudly. Each error names the offending
  rule id and gives per-type guidance: Sat → "use MM instead"
  (mirroring NFsim's own policy at `NFinput.cpp:2459`); Hill →
  "use generate_network + ODE"; FunctionProduct → "rewrite as a
  single Function".

- **Feature-coverage tests for MM and TFUN** (`6421017`,
  `cb03071`).  `ft_mm_ratelaw.bngl` exercises `MM(kcat,Km)`
  against NFsim; `ft_tfun.bngl` exercises the new lowercase
  `tfun()` syntax against BNG dev-branch ODE (set
  `BNG2=$HOME/Code/bionetgen/bng2/BNG2.pl` to regenerate). Both
  code paths now have regression coverage; both exposed real RM
  bugs while being authored.

- **`examples/embed.cpp` + `examples/CMakeLists.txt`**
  (`a797b81`).  Minimal compilable C++ embedding example showing
  both stateless `run()` and stateful session usage. Off by
  default; opt in via `cmake -DRULEMONKEY_BUILD_EXAMPLES=ON`.

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

- **README**: replaced the bare "Public API" snippet with a longer
  "Embedding (C++ API)" section that pairs stateless and session
  usage and points at the new example and public headers
  (`a797b81`).

- **`harness/benchmark_feature_coverage.py`**: added
  `_copy_aux_files()` to stage `*.tfun` next to the BNGL when
  invoking BNG2 / NFsim from a tempdir; `_run_one_nfsim_rep` now
  runs with `cwd=tmpdir` so NFsim's relative-path lookups
  resolve. New `RM_ONLY` set (currently empty) for any future
  model that has no third-party reference.

- **Two existing models cleaned up post-`edg_*` benchmark**
  (`6617a84`): `edg_fixed_competition` dropped a bogus
  `conserved Total_S = 90` invariant (S is consumed by catalysis,
  not conserved); `edg_ring_break_constraint` corrected
  `conserved A_total` from 30 → 40 (miscounted seed: 10×2-mer +
  5×4-mer = 40 A's); `edg_deep_param_chain` magnitudes bumped so
  the default 5-rep ODE comparison is stable. Harness routes
  `edg_time_dependent_rate` and `edg_deep_param_chain` through
  ODE verdict (NFsim rejects `time()` and function-of-function
  chains).

- `tests/reference/basicmodels/PROVENANCE.md` rewritten to lead
  with what the suite *is* (29 imported tests, source, reference
  flow) and treat the seven upstream NFsim tests not carried over
  as a clearly-labeled appendix grouped by reason.
  `harness/basicmodels.py` and
  `harness/generate_basicmodels_refs.py` docstrings tightened to
  a one-line origin pointer.

### Removed

- **Dead `extra_eff` machinery** (`d061f18`). After the self-
  binding propensity fix above, `compute_extra()`,
  `PerMolRuleData::extra`, `RuleState::total_extra`, and
  `extra_eff` in `compute_propensity` became unused — within-mol
  pair removal is now handled at sample time by
  `select_reactants`, not at propensity time. 35 lines removed
  from `engine.cpp`.

- **Four upstream NFsim regression tests dropped from the
  basicmodels suite** (`85feae1`, `9fb2efb`). All four tested
  NFsim-specific behaviors that don't apply to RuleMonkey:
  - `r33` and `r35` pinned NFsim issues #22 / #21 ("occupied-
    site bond error") and #14 (RHS `.` between products that
    NFsim splits anyway). On the BNGL these tests carry,
    BNG2.pl's `generate_network` produces the chemistry-correct
    behavior (bound by free-B count for r33; zero reactions for
    r35) and RuleMonkey matches BNG2.pl. The NFsim references
    captured the historic NFsim quirks, which by design diverge
    from BNGL strict.
  - `r31` is a crash regression test with no `begin observables`
    block (the author's own comment: *"validation harness will
    run NFsim on this XML and ensure it doesn't crash"*).
  - `r34` includes a `begin observables` block but the author
    deliberately commented out the only line in it.

  After the removals the suite is clean **29 PASS / 0 FAIL /
  0 NO_MATCH**. Joins the pre-existing r27 / r28 / r36 in the
  PROVENANCE appendix.

### Verification

End-to-end benchmark state on 2026-04-29 (post-fix):

| Suite                                       | Result |
|---------------------------------------------|---|
| `feature_coverage` (77 models, --reps 5)    | 77 PASS / 0 FAIL |
| `benchmark_full --tier full` (71 corpus)    | 71 PASS / 0 FAIL |
| `nfsim_basicmodels`         (29 models)     | 29 PASS / 0 FAIL |

177 / 177 models PASS RM-vs-NFsim z-score (or RM-vs-ODE rel-err
for `NFSIM_UNRELIABLE` models) post-fix.

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
  - `feature_coverage/` — 52 BNGL feature-coverage models with invariants
    and golden values.
  - `corpus/` — 71 real-world rule-based models for efficiency and
    correctness benchmarking.
  - `nfsim_basicmodels/` — 36 models from the NFsim parity suite.
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

[3.1.0]: https://github.com/wshlavacek/RuleMonkey/releases/tag/v3.1.0
[3.0.0]: https://github.com/wshlavacek/RuleMonkey/releases/tag/v3.0.0
