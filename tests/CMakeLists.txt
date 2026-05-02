# Smoke test: parses an XML model and runs a short simulation.
add_executable(smoke_test cpp/smoke_test.cpp)
target_link_libraries(smoke_test PRIVATE rulemonkey)

add_test(
    NAME smoke_test
    COMMAND smoke_test ${CMAKE_SOURCE_DIR}/tests/reference/nfsim/xml/A_plus_A.xml
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
        ${CMAKE_SOURCE_DIR}/tests/reference/nfsim/xml/A_plus_A.xml
        ${CMAKE_SOURCE_DIR}/tests/models/feature_coverage/xml/ft_mm_ratelaw.xml
        ${CMAKE_SOURCE_DIR}/tests/cpp/derived_param_model.xml
        ${CMAKE_SOURCE_DIR}/tests/cpp/out_of_order_param_model.xml
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
        ${CMAKE_SOURCE_DIR}/tests/reference/nfsim/xml/A_plus_A.xml
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
        ${CMAKE_SOURCE_DIR}/tests/reference/nfsim/xml/A_plus_A.xml
        ${CMAKE_SOURCE_DIR}/tests/models/feature_coverage/xml/ft_mm_ratelaw.xml
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
        ${CMAKE_SOURCE_DIR}/tests/reference/nfsim/xml/A_plus_A.xml
)

# Direct unit tests for the BNGL math expression parser/evaluator.  Pin
# down review #1 (unary minus precedence), #2 (&& / || short-circuit),
# and #6 (== / != bit-exact) at the parse-and-evaluate level so a
# regression surfaces here as a one-line miscompare instead of a
# trajectory drift on some incidental BNGL parity model.  Needs PRIVATE
# access to cpp/rulemonkey/expr_eval.hpp (not part of the public API).
add_executable(expr_eval_test cpp/expr_eval_test.cpp)
target_link_libraries(expr_eval_test PRIVATE rulemonkey)
target_include_directories(expr_eval_test PRIVATE ${CMAKE_SOURCE_DIR}/cpp/rulemonkey)

add_test(NAME expr_eval_test COMMAND expr_eval_test)

# Direct unit tests for TableFunction (TFUN) interpolation.  Pin down
# review #7 (step is right-continuous, matches BNG Network3) along with
# linear interp + boundary clamp + ctor validation.  Same private-header
# access pattern as expr_eval_test.
add_executable(table_function_test cpp/table_function_test.cpp)
target_link_libraries(table_function_test PRIVATE rulemonkey)
target_include_directories(table_function_test PRIVATE ${CMAKE_SOURCE_DIR}/cpp/rulemonkey)

add_test(NAME table_function_test COMMAND table_function_test)
