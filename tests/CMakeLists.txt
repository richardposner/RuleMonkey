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
