# Smoke test: parses an XML model and runs a short simulation.
add_executable(smoke_test cpp/smoke_test.cpp)
target_link_libraries(smoke_test PRIVATE rulemonkey)

add_test(
    NAME smoke_test
    COMMAND smoke_test ${CMAKE_SOURCE_DIR}/tests/reference/nfsim/xml/A_plus_A.xml
)
