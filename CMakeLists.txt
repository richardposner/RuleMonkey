cmake_minimum_required(VERSION 3.25)

project(RuleMonkey VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

if(APPLE)
    set(CMAKE_OSX_ARCHITECTURES "arm64;x86_64" CACHE STRING "Universal binary")
endif()

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
endif()

# ---------------------------------------------------------------------------
# librulemonkey: the core simulation engine
# ---------------------------------------------------------------------------
add_library(rulemonkey STATIC
    cpp/rulemonkey/engine.cpp
    cpp/rulemonkey/simulator.cpp
    cpp/rulemonkey/expr_eval.cpp
    cpp/rulemonkey/table_function.cpp
)

target_include_directories(rulemonkey
    PUBLIC  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/rulemonkey
)

target_compile_features(rulemonkey PUBLIC cxx_std_17)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(rulemonkey PRIVATE
        -Wall -Wextra -Wno-deprecated-declarations
    )
endif()

# ---------------------------------------------------------------------------
# rm_driver: batch CLI
# ---------------------------------------------------------------------------
add_executable(rm_driver cpp/cli/rm_driver.cpp)
target_link_libraries(rm_driver PRIVATE rulemonkey)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(rm_driver PRIVATE -Wno-deprecated-declarations)
endif()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
enable_testing()
add_subdirectory(tests)
