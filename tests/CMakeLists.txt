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

# Direct unit tests for the BNGL math expression parser/evaluator.  Pin
# down review #1 (unary minus precedence), #2 (&& / || short-circuit),
# and #6 (== / != bit-exact) at the parse-and-evaluate level so a
# regression surfaces here as a one-line miscompare instead of a
# trajectory drift on some incidental BNGL parity model.  Needs PRIVATE
# access to cpp/rulemonkey/expr_eval.hpp (not part of the public API).
add_executable(expr_eval_test cpp/expr_eval_test.cpp)
target_link_libraries(expr_eval_test PRIVATE rulemonkey)
target_include_directories(expr_eval_test PRIVATE ${_rulemonkey_private_include_dir})

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
