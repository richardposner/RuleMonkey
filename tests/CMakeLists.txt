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
)
