# Smoke test: parses an XML model and runs a short simulation.
set(_rulemonkey_test_dir "${PROJECT_SOURCE_DIR}/tests")
set(_rulemonkey_private_include_dir "${PROJECT_SOURCE_DIR}/cpp/rulemonkey")

add_executable(smoke_test cpp/smoke_test.cpp)
target_link_libraries(smoke_test PRIVATE rulemonkey)

add_test(
    NAME smoke_test
    COMMAND smoke_test ${_rulemonkey_test_dir}/reference/nfsim/xml/A_plus_A.xml
)

# set_param regression: verifies overrides reach the engine for Ele
# rate constants, MM rate constants, and initial species
# concentrations.  Without the apply_overrides re-resolve fix these
# paths silently ignored set_param.
add_executable(set_param_test cpp/set_param_test.cpp)
target_link_libraries(set_param_test PRIVATE rulemonkey)

add_test(
    NAME set_param_test
    COMMAND set_param_test
        ${_rulemonkey_test_dir}/reference/nfsim/xml/A_plus_A.xml
        ${_rulemonkey_test_dir}/models/feature_coverage/xml/ft_mm_ratelaw.xml
        ${_rulemonkey_test_dir}/cpp/derived_param_model.xml
        ${_rulemonkey_test_dir}/cpp/out_of_order_param_model.xml
        ${_rulemonkey_test_dir}/models/feature_coverage/xml/ft_energy_arrhenius.xml
        ${_rulemonkey_test_dir}/reference/nfsim/xml/bench_blbr_cooperativity_posner2004.xml
)

# Corpus-wide guard (issue #23): a no-op set_param switches the parameter
# cascade onto the `<Parameter expr=>` source, and no model's parameters
# may drift off their loaded `<Parameter value=>` when it does.  Load-only,
# so the whole 180+ model tree sweeps in about a second.
add_executable(param_expr_gate_test cpp/param_expr_gate_test.cpp)
target_link_libraries(param_expr_gate_test PRIVATE rulemonkey)

add_test(
    NAME param_expr_gate_test
    COMMAND param_expr_gate_test ${_rulemonkey_test_dir}
)

# Seed-species amount API (issue #23): introspection plus pin/release of
# the initial population, so a driver can walk a dose-response curve on
# one loaded model instead of re-baking an XML per point.  Covers the two
# shapes set_param cannot reach — a literal <Species concentration=> with
# no parameter behind it, and a Fixed="1" clamp target.
add_executable(initial_amount_test cpp/initial_amount_test.cpp)
target_link_libraries(initial_amount_test PRIVATE rulemonkey)

add_test(
    NAME initial_amount_test
    COMMAND initial_amount_test
        ${_rulemonkey_test_dir}/reference/nfsim/xml/A_plus_A.xml
        ${_rulemonkey_test_dir}/models/feature_coverage/xml/ft_clamped_species.xml
)

# Homodimer rate parity: small-N steady state of A+A<->AA must match the
# chemical-master-equation analytic mean.  Pins down the same_components
# branch in compute_propensity (deflated `(ab²−ab_sq)/2`) plus the
# retry-until-distinct sampler against the CME — a regression in either
# would shift the equilibrium population by O(1/N).
add_executable(homodimer_rate_test cpp/homodimer_rate_test.cpp)
target_link_libraries(homodimer_rate_test PRIVATE rulemonkey)

add_test(
    NAME homodimer_rate_test
    COMMAND homodimer_rate_test
        ${_rulemonkey_test_dir}/reference/nfsim/xml/A_plus_A.xml
)

# DOR1 rate parity (issue #34): a bimolecular rule with ONE tagged reactant
# and a local-function rate.  No branch of the propensity code covered that
# shape, so the rule fired once and then read a never-populated accumulator
# as its propensity — silently inert.  Four arms (tag on either reactant,
# molecule vs complex observable scope, single- vs multi-molecule tagged
# pattern) are checked against an exact binomial death-process oracle.
add_executable(local_fcn_bimol_test cpp/local_fcn_bimol_test.cpp)
target_link_libraries(local_fcn_bimol_test PRIVATE rulemonkey)

add_test(
    NAME local_fcn_bimol_test
    COMMAND local_fcn_bimol_test
        ${_rulemonkey_test_dir}/models/feature_coverage/xml/ft_local_fcn_bimol.xml
)

# symmetry_factor on an MM(kcat,Km) rate law (issue #37): every rate law but
# the Michaelis-Menten branch applied BNG2's factor, and no corpus model
# paired one with MM.  The factor scales the substrate COUNT inside the law,
# not the finished propensity — MM is nonlinear in that count, so the two
# differ everywhere except below saturation.  Three arms (symmetric dimer
# saturated, symmetric dimer linear, asymmetric control) against a mean-field
# firing-count oracle; no single regime separates all three candidates.
add_executable(mm_symmetry_test cpp/mm_symmetry_test.cpp)
target_link_libraries(mm_symmetry_test PRIVATE rulemonkey)

add_test(
    NAME mm_symmetry_test
    COMMAND mm_symmetry_test
        ${_rulemonkey_test_dir}/cpp/mm_symmetry_model.xml
)

# The Km <= 0 boundary of MM(kcat,Km) (issue #46).  The division guard
# answered 0 at Km == 0, where the limit is kcat*min(S,E), so the rule went
# silently inert (and a model starting at S > E froze the moment its
# substrate decayed to S == E); the naive sFree cancelled catastrophically on
# the approach (17% fast at Km=1e-14, dead below ~1e-15); and a negative Km
# passed its NaN through as the propensity, which the sampler then drew
# events against.  Closed-form oracles, no external simulator.
add_executable(mm_km_domain_test cpp/mm_km_domain_test.cpp)
target_link_libraries(mm_km_domain_test PRIVATE rulemonkey)

add_test(
    NAME mm_km_domain_test
    COMMAND mm_km_domain_test
        ${_rulemonkey_test_dir}/models/feature_coverage/xml/ft_mm_ratelaw.xml
)

# Load-time diagnostics for the two MM(kcat,Km) constructs RM reads one way
# and cannot verify (issue #45): a substrate pattern matching more than one
# species, where BNG2's network expansion gives each matching species its own
# copy of the enzyme pool, and a symmetry_factor on a rule that transforms
# both reactant patterns, where the scalar cannot be attributed to a slot.
# Both are Warn severity; mm_symmetry_model is the negative control.
add_executable(mm_diagnostics_test cpp/mm_diagnostics_test.cpp)
target_link_libraries(mm_diagnostics_test PRIVATE rulemonkey)

add_test(
    NAME mm_diagnostics_test
    COMMAND mm_diagnostics_test
        ${_rulemonkey_test_dir}/cpp/mm_diagnostics_model.xml
        ${_rulemonkey_test_dir}/cpp/mm_symmetry_model.xml
)

# Observable scope inside a local function (issue #38): an observable applied
# to the function's tag is local, a bare one is global, and BNG2's XML marks
# both identically — the split only survives in the <Expression>.  RM
# localized every one of them, so a bare observable read 0 and its rule went
# silently inert.  Covers the unimolecular and bimolecular propensity paths,
# one observable used at both scopes across two functions (which the union of
# per-function local sets cannot express), a bare observable whose value moves
# mid-run, and the warning for both scopes inside a single function.
add_executable(local_fcn_global_obs_test cpp/local_fcn_global_obs_test.cpp)
target_link_libraries(local_fcn_global_obs_test PRIVATE rulemonkey)

add_test(
    NAME local_fcn_global_obs_test
    COMMAND local_fcn_global_obs_test
        ${_rulemonkey_test_dir}/models/feature_coverage/xml/ft_local_fcn_global_obs.xml
        ${_rulemonkey_test_dir}/models/feature_coverage/xml/ft_local_fcn_mixed_scope.xml
        ${_rulemonkey_test_dir}/cpp/local_fcn_mixed_ref_model.xml
)

# `reactant_N()` in a rate law (issue #59): the empty-bodied placeholder a
# BNGL model declares when it wants the match count of a rule's Nth reactant
# pattern inside that rule's own rate law.  RM read it as the constant zero,
# so any rate law built on one evaluated to zero and the rule never fired —
# no exception, no warning, just a run with no events.  Pins the reported
# rate value, the propensity behind it (measured as the firing fraction at
# the mean waiting time), the coverage model's placeholder arms against
# their observable-only controls, and the two shapes RM refuses at load.
add_executable(reactant_count_rate_test cpp/reactant_count_rate_test.cpp)
target_link_libraries(reactant_count_rate_test PRIVATE rulemonkey)

add_test(
    NAME reactant_count_rate_test
    COMMAND reactant_count_rate_test
        ${_rulemonkey_test_dir}/cpp/reactant_count_nucleation_model.xml
        ${_rulemonkey_test_dir}/models/feature_coverage/xml/ft_reactant_count_rate.xml
        ${_rulemonkey_test_dir}/cpp/reactant_count_arity_model.xml
        ${_rulemonkey_test_dir}/cpp/reactant_count_local_fcn_model.xml
        ${_rulemonkey_test_dir}/cpp/reactant_count_total_rate_model.xml
)

# Cost of that rescan (issue #40): #38 made a local rate reading ANY global
# rescan every instance after every event, when what gates an O(N) rescan is
# whether the global actually moved — and a bare observable is usually a
# volume, a total or a clock that never does.  Runs a 22000-instance model
# whose bare observable is constant against a twin with that observable
# folded into a parameter, and requires the two to agree on the trajectory
# (the skip is exact) AND on the wall time (the skip is free).  The twin also
# covers the other direction: add_molecules moves the observable from outside
# the event loop, where the value cache is the only thing watching.
add_executable(local_fcn_rescan_gate_test cpp/local_fcn_rescan_gate_test.cpp)
target_link_libraries(local_fcn_rescan_gate_test PRIVATE rulemonkey)

add_test(
    NAME local_fcn_rescan_gate_test
    COMMAND local_fcn_rescan_gate_test
        ${_rulemonkey_test_dir}/models/feature_coverage/xml/ss_local_fcn_const_global.xml
        ${_rulemonkey_test_dir}/cpp/local_fcn_const_global_folded_model.xml
)

# save_state / load_state regression: round-trip a mid-run snapshot
# back into a fresh simulator and verify (a) the continuation matches
# an uninterrupted run, (b) loading state from a structurally
# different XML is refused via the schema fingerprint, and (c) the
# pre-V2 marker is refused outright.
add_executable(save_load_test cpp/save_load_test.cpp)
target_link_libraries(save_load_test PRIVATE rulemonkey)

add_test(
    NAME save_load_test
    COMMAND save_load_test
        ${_rulemonkey_test_dir}/reference/nfsim/xml/A_plus_A.xml
        ${_rulemonkey_test_dir}/models/feature_coverage/xml/ft_mm_ratelaw.xml
)

# Public-API coverage: surface-area regressions on getters and behavioral
# methods that aren't otherwise exercised by the parity tests.  Pins down
# step_to / has_session / get_observable_values / parameter_names /
# observable_names / xml_path / method / unsupported_features and the
# behavioral effect of set_block_same_complex_binding / set_molecule_limit.
add_executable(api_coverage_test cpp/api_coverage_test.cpp)
target_link_libraries(api_coverage_test PRIVATE rulemonkey)

add_test(
    NAME api_coverage_test
    COMMAND api_coverage_test
        ${_rulemonkey_test_dir}/reference/nfsim/xml/A_plus_A.xml
)

# Direct unit tests for RuleMonkey's expression layer (issue #6): the
# expr::collect_variables dependency scanner and the ExprTk semantics RM
# depends on (unary-minus precedence, ^ associativity, bit-exact == / !=,
# lazy if).  A regression surfaces here as a one-line miscompare instead
# of a trajectory drift on some incidental BNGL parity model.  Needs
# PRIVATE access to cpp/rulemonkey/expr_eval.hpp (not part of the public
# API) and to the vendored bngsim/expression.hpp evaluator header.
add_executable(expr_eval_test cpp/expr_eval_test.cpp)
target_link_libraries(expr_eval_test PRIVATE rulemonkey)
target_include_directories(expr_eval_test PRIVATE
    ${_rulemonkey_private_include_dir}
    ${PROJECT_SOURCE_DIR}/third_party/bngsim_expr/include)

add_test(NAME expr_eval_test COMMAND expr_eval_test)

# Direct unit tests for TableFunction (TFUN) interpolation.  Pin down
# review #7 (step is right-continuous, matches BNG Network3) along with
# linear interp + boundary clamp + ctor validation.  Same private-header
# access pattern as expr_eval_test.
add_executable(table_function_test cpp/table_function_test.cpp)
target_link_libraries(table_function_test PRIVATE rulemonkey)
target_include_directories(table_function_test PRIVATE ${_rulemonkey_private_include_dir})

add_test(NAME table_function_test COMMAND table_function_test)

# Direct unit tests for the canonical complex labeler (issue #9 §2;
# dev/canonical_labeling_plan_2026_05_16.md steps 1-2).  Hand-built
# ComplexGraph shapes pinned to exact canonical strings, plus a
# property-based test that applies random graph isomorphisms (molecule
# + interchangeable-component permutations) and asserts the canonical
# label is invariant — the test that catches refinement / render /
# bond-labeling bugs.  Needs PRIVATE access to cpp/rulemonkey/canonical.hpp
# (the labeler is engine-internal, not part of the public API).
add_executable(canonical_test cpp/canonical_test.cpp)
target_link_libraries(canonical_test PRIVATE rulemonkey)
target_include_directories(canonical_test PRIVATE ${_rulemonkey_private_include_dir})

add_test(NAME canonical_test COMMAND canonical_test)

# eBNGL energy-rule expansion: pins the Sekar port (energy_expand.{hpp,cpp})
# to the reference NFsim's pre-computed rate constants for the v40 test-suite
# model, plus the cooperative-scaffold always/conditional classification.
add_executable(energy_expand_test cpp/energy_expand_test.cpp)
target_link_libraries(energy_expand_test PRIVATE rulemonkey)
target_include_directories(energy_expand_test PRIVATE ${_rulemonkey_private_include_dir})

add_test(NAME energy_expand_test COMMAND energy_expand_test)

# Runtime pattern-parser unit test (issue #9 §1 step 2).  Pure unit
# test: parse_species_pattern is string + Model -> Pattern with no
# engine state, so a hand-built Model drives the accepted grammar,
# model resolution (types / components / states / symmetric components
# / bonds) and the exact-species error surface.  Needs PRIVATE access
# to cpp/rulemonkey/pattern_parser.hpp (the parser is engine-internal,
# not part of the public API).
add_executable(pattern_parser_test cpp/pattern_parser_test.cpp)
target_link_libraries(pattern_parser_test PRIVATE rulemonkey)
target_include_directories(pattern_parser_test PRIVATE ${_rulemonkey_private_include_dir})

add_test(NAME pattern_parser_test COMMAND pattern_parser_test)

# Pattern-keyed species-method integration test (issue #9 §1 step 3):
# drives get_species_count / add_species / remove_species /
# set_species_count through the public API — counts agree with
# enumerate_species(), single- and multi-molecule add/remove round-trip,
# the SSA loop still runs after a structural mutation, and the error
# surface (no session, malformed pattern, bad counts) throws.
add_executable(species_methods_test cpp/species_methods_test.cpp)
target_link_libraries(species_methods_test PRIVATE rulemonkey)

add_test(
    NAME species_methods_test
    COMMAND species_methods_test
        ${_rulemonkey_test_dir}/reference/nfsim/xml/A_plus_A.xml
        ${_rulemonkey_test_dir}/cpp/species_order_model.xml
)

# Species enumeration + `.species` output (issue #9 §2): exercises the
# batch pipeline end to end — pool walk, extract_complex, canonical
# labeling, isomorphism dedup, sorted SpeciesRow rows, and the
# BNG-format `.species` file writer — through the public simulator API.
# Also drives the cached-incremental label layer (plan §6 step 5): in
# Debug/ASan builds enumerate_species cross-checks every cached label
# against a from-scratch recompute, so the merge/split/ring models swept
# repeatedly here abort on a missed cache-invalidation hook.
add_executable(species_enumeration_test cpp/species_enumeration_test.cpp)
target_link_libraries(species_enumeration_test PRIVATE rulemonkey)

add_test(
    NAME species_enumeration_test
    COMMAND species_enumeration_test
        ${_rulemonkey_test_dir}/reference/nfsim/xml/A_plus_A.xml
        ${_rulemonkey_test_dir}/models/feature_coverage/xml/ss_symmetric_homopoly.xml
        ${_rulemonkey_test_dir}/models/feature_coverage/xml/ft_ring_closure.xml
)

# Error-path regression: pins down that the documented public-API
# error surfaces (missing file, malformed XML, unknown set_param,
# mutators called while a session is active) throw std::runtime_error
# rather than silently failing or crashing.  Uses a tmp file path
# under the build dir for the malformed-XML case.
add_executable(error_paths_test cpp/error_paths_test.cpp)
target_link_libraries(error_paths_test PRIVATE rulemonkey)

add_test(
    NAME error_paths_test
    COMMAND error_paths_test
        ${_rulemonkey_test_dir}/reference/nfsim/xml/A_plus_A.xml
        ${CMAKE_CURRENT_BINARY_DIR}/error_paths_malformed.xml
)

# Negative-propensity clamp diagnostic regression: pins down that
# function rate laws evaluating to a negative number (a) get clamped
# to zero rather than aborting the run, and (b) emit exactly one WARN
# line per rule on first occurrence — the diagnostic surface NFsim
# refuses the model entirely on.  Without the per-rule latch a noisy
# oscillator would spam the log; without the WARN at all the author
# would have no signal that the BNGL is mis-authored for SSA.
add_executable(negative_rate_clamp_test cpp/negative_rate_clamp_test.cpp)
target_link_libraries(negative_rate_clamp_test PRIVATE rulemonkey)

add_test(
    NAME negative_rate_clamp_test
    COMMAND negative_rate_clamp_test
        ${_rulemonkey_test_dir}/cpp/negative_rate_clamp_model.xml
        ${CMAKE_CURRENT_BINARY_DIR}/negative_rate_clamp_stderr.txt
)

# N-ary reactant rules (issues #24, #26): a rule such as `A + A + A -> P` used
# to collapse its 2nd and 3rd patterns into the engine's single slot-B tail,
# score zero embeddings, and never fire — with mass still conserved and
# nothing reported, so the trajectory looked valid.  Covers that the
# reproducer now fires and balances mass, that the realized rate matches the
# closed-form k·sf·N(N-1)(N-2) on a constant-population fixture (the check a
# wrong symmetry factor fails), that a `.`-joined reactant complex is
# accepted and its overlapping matches neither double-consume a molecule nor
# skew the rate (#26), and that shapes outside the implemented path are still
# refused rather than going quietly inert again.
add_executable(trimolecular_test cpp/trimolecular_test.cpp)
target_link_libraries(trimolecular_test PRIVATE rulemonkey)

add_test(
    NAME trimolecular_test
    COMMAND trimolecular_test
        ${_rulemonkey_test_dir}/cpp/trimolecular_model.xml
        ${_rulemonkey_test_dir}/cpp/trimolecular_rate_model.xml
        ${_rulemonkey_test_dir}/cpp/trimolecular_multimol_model.xml
        ${_rulemonkey_test_dir}/cpp/trimolecular_overlap_model.xml
        ${_rulemonkey_test_dir}/cpp/trimolecular_deadend_model.xml
        ${_rulemonkey_test_dir}/cpp/trimolecular_overlap_delete_model.xml
        ${_rulemonkey_test_dir}/cpp/trimolecular_disjoint_model.xml
        ${_rulemonkey_test_dir}/cpp/trimolecular_toomany_model.xml
        ${_rulemonkey_test_dir}/cpp/trimolecular_ratelaw_model.xml
)

# Bimolecular multi-molecule reactant rate parity: the sampler weights the
# seed molecule by the number of seed embeddings that extend to a whole match
# (count_multi_mol_fast's c) but used to choose among all S >= c of them,
# turning the dead ends into null events and firing the rule at c/S of its
# rate.  The fixture makes c/S = 1/2 — an A(d,d) bonded to a D and an E, with
# the pattern A(d!1).D(d!1) — and holds every population constant so the
# closed-form propensity applies for the whole run.  The second fixture
# covers the other multi-molecule hazard in the same sampler: a molecule
# reached through the pattern's bonds also matching the other slot, which
# under DeleteMolecules deleted it twice and broke mass balance.
add_executable(multimol_rate_test cpp/multimol_rate_test.cpp)
target_link_libraries(multimol_rate_test PRIVATE rulemonkey)

add_test(
    NAME multimol_rate_test
    COMMAND multimol_rate_test
        ${_rulemonkey_test_dir}/cpp/multimol_deadend_model.xml
        ${_rulemonkey_test_dir}/cpp/multimol_overlap_model.xml
)

# Global-function exposure regression (issue #7): pins down that
# Result::function_names / function_data, RuleMonkeySimulator::function_names()
# and get_function_values() expose BNGL global-function values, that the
# per-sample columns stay algebraically consistent with the observables
# they derive from (including nested function-of-function settle order),
# and that a model with no functions presents an empty — not absent —
# function surface.
add_executable(function_values_test cpp/function_values_test.cpp)
target_link_libraries(function_values_test PRIVATE rulemonkey)

add_test(
    NAME function_values_test
    COMMAND function_values_test
        ${_rulemonkey_test_dir}/models/feature_coverage/xml/ft_nested_functions.xml
        ${_rulemonkey_test_dir}/reference/nfsim/xml/A_plus_A.xml
)

# evaluate_expression session API (issue #9 §1): pins down that an
# arbitrary BNGL expression evaluated against a live session resolves
# parameters, observables, global functions, and the time() / t clock
# to the same values the dedicated accessors report, that the `extra`
# map shadows model symbols without mutating session state, and that a
# compile failure or a no-session call throws.
add_executable(evaluate_expression_test cpp/evaluate_expression_test.cpp)
target_link_libraries(evaluate_expression_test PRIVATE rulemonkey)

add_test(
    NAME evaluate_expression_test
    COMMAND evaluate_expression_test
        ${_rulemonkey_test_dir}/models/feature_coverage/xml/ft_nested_functions.xml
)

# rm_driver `.gdat` function-column opt-in (issue #7): the global-function
# columns must appear ONLY with --print-functions, matching BNGL's
# `print_functions=>1`.  Two CTest fixtures on the same function-bearing
# model assert the header content directly via PASS/FAIL regexes — no C++
# harness needed.  Guarded on the rm_driver target so a tests-only build
# (RULEMONKEY_BUILD_CLI=OFF) still configures.
if(TARGET rm_driver)
    add_test(
        NAME rm_driver_functions_default_off
        COMMAND rm_driver
            ${_rulemonkey_test_dir}/models/feature_coverage/xml/ft_nested_functions.xml
            10 5 7
    )
    # Default run: observables present, function columns absent.
    set_tests_properties(rm_driver_functions_default_off PROPERTIES
        PASS_REGULAR_EXPRESSION "S_total"
        FAIL_REGULAR_EXPRESSION "frac_P")

    add_test(
        NAME rm_driver_functions_print_flag
        COMMAND rm_driver
            ${_rulemonkey_test_dir}/models/feature_coverage/xml/ft_nested_functions.xml
            10 5 7 --print-functions
    )
    # With --print-functions: all three global functions appear in the header.
    set_tests_properties(rm_driver_functions_print_flag PROPERTIES
        PASS_REGULAR_EXPRESSION "frac_P.*modulator.*phos_rate")
endif()

# Cooperative-cancellation regression: pins down the CancelCallback
# contract on run() / simulate() / step_to() — pre-cancelled callback
# throws immediately, Cancelled inherits std::runtime_error, mid-session
# cancellation leaves the session at a partial time and is recoverable,
# and an always-true callback produces a bit-identical trajectory.
# Embedders (notably BNGsim's `timeout` kwarg path) rely on all four.
add_executable(cancellation_test cpp/cancellation_test.cpp)
target_link_libraries(cancellation_test PRIVATE rulemonkey)

add_test(
    NAME cancellation_test
    COMMAND cancellation_test
        ${_rulemonkey_test_dir}/reference/nfsim/xml/A_plus_A.xml
)

# parameter_scan / bifurcate regression (issue #8): pins down value
# generation (explicit list, linear + geometric ranges, degenerate
# single-point / zero-width cases), endpoint shape, fixed-seed
# determinism, the monotone AA response of A+A<->AA to the forward
# rate, global-function columns, the ascending-axis alignment of
# bifurcate's forward/backward branches, ScanSpec validation, and that
# a sweep restores the swept parameter's override state.
add_executable(parameter_scan_test cpp/parameter_scan_test.cpp)
target_link_libraries(parameter_scan_test PRIVATE rulemonkey)

add_test(
    NAME parameter_scan_test
    COMMAND parameter_scan_test
        ${_rulemonkey_test_dir}/reference/nfsim/xml/A_plus_A.xml
        ${_rulemonkey_test_dir}/models/feature_coverage/xml/ft_nested_functions.xml
)

# Explicit output-time (sample_times) regression (issue #16): pins down
# that TimeSpec::sample_times overrides the uniform n_points grid — an
# explicit subset of the uniform grid is bit-identical to the full run at
# the shared instants (sampling is non-invasive and event_count is
# unchanged), arbitrary non-uniform instants are echoed verbatim with one
# row each, and an out-of-order list throws.  Uses A+A<->AA so the
# observable values are exact integer counts and the bit-identical check is
# a true equality test.
add_executable(sample_times_test cpp/sample_times_test.cpp)
target_link_libraries(sample_times_test PRIVATE rulemonkey)

add_test(
    NAME sample_times_test
    COMMAND sample_times_test
        ${_rulemonkey_test_dir}/reference/nfsim/xml/A_plus_A.xml
)

# Session counterpart of sample_times_test: the stateful simulate(const
# TimeSpec&) overload honors TimeSpec::sample_times bit-identically to a
# uniform session run at the shared instants, works mid-protocol (continuing
# the trajectory after a prior segment), and validates that the segment start
# matches the current session time.
add_executable(session_sample_times_test cpp/session_sample_times_test.cpp)
target_link_libraries(session_sample_times_test PRIVATE rulemonkey)

add_test(
    NAME session_sample_times_test
    COMMAND session_sample_times_test
        ${_rulemonkey_test_dir}/reference/nfsim/xml/A_plus_A.xml
)

# rm_scan CLI smoke + `.scan` header regressions (issue #8).  Like the
# rm_driver function-column tests, these assert output shape directly
# via PASS/FAIL regexes.  Guarded on the rm_scan target so a tests-only
# build (RULEMONKEY_BUILD_CLI=OFF) still configures.
if(TARGET rm_scan)
    # parameter_scan: header is `# <param>` + observable names; no
    # function columns without --print-functions.
    add_test(
        NAME rm_scan_parameter_scan
        COMMAND rm_scan
            ${_rulemonkey_test_dir}/models/feature_coverage/xml/ft_nested_functions.xml
            --param k0 --min 0.01 --max 0.2 --n-points 4 --t-end 20
    )
    set_tests_properties(rm_scan_parameter_scan PROPERTIES
        PASS_REGULAR_EXPRESSION "# k0\t.*S_total"
        FAIL_REGULAR_EXPRESSION "phos_rate")

    # --print-functions appends the global-function columns.
    add_test(
        NAME rm_scan_print_functions
        COMMAND rm_scan
            ${_rulemonkey_test_dir}/models/feature_coverage/xml/ft_nested_functions.xml
            --param k0 --min 0.01 --max 0.2 --n-points 4 --t-end 20 --print-functions
    )
    set_tests_properties(rm_scan_print_functions PROPERTIES
        PASS_REGULAR_EXPRESSION "frac_P.*modulator.*phos_rate")

    # --bifurcate: every observable gets paired `_fwd` / `_bwd` columns.
    add_test(
        NAME rm_scan_bifurcate
        COMMAND rm_scan
            ${_rulemonkey_test_dir}/reference/nfsim/xml/A_plus_A.xml
            --param kp --min 0.001 --max 0.006 --n-points 5 --t-end 20 --bifurcate
    )
    set_tests_properties(rm_scan_bifurcate PROPERTIES
        PASS_REGULAR_EXPRESSION "A_1_fwd\tA_1_bwd")

    # An invalid spec (n_points == 1 with distinct bounds) exits non-zero.
    add_test(
        NAME rm_scan_rejects_bad_spec
        COMMAND rm_scan
            ${_rulemonkey_test_dir}/reference/nfsim/xml/A_plus_A.xml
            --param kp --min 0.001 --max 0.006 --n-points 1 --t-end 20
    )
    set_tests_properties(rm_scan_rejects_bad_spec PROPERTIES WILL_FAIL TRUE)
endif()

# BNG2-oracle suite (tests/bng_oracle): models whose reference is BioNetGen's
# own network expansion rather than an NFsim ensemble, because they pin
# behaviour NFsim does not have.  Scored by z against the committed BNG2 mean;
# see tests/bng_oracle/README.md for why the suite is separate and how the two
# reference kinds (ode / ssa) are chosen.
#
# Needs rm_driver and a python3; skipped when either is absent so a
# library-only or CLI-disabled configure still passes.
if(TARGET rm_driver)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
        add_test(
            NAME bng_oracle
            COMMAND ${Python3_EXECUTABLE} ${_rulemonkey_test_dir}/bng_oracle/check.py
                --driver $<TARGET_FILE:rm_driver>
        )
        # The 400-seed sampler model dominates; comfortably inside this on a
        # laptop, headroom for a loaded CI runner.
        set_tests_properties(bng_oracle PROPERTIES TIMEOUT 600)
    endif()
endif()
