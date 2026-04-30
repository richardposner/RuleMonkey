cmake_minimum_required(VERSION 3.25)

project(RuleMonkey VERSION 3.1.0 LANGUAGES CXX)

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
            $<INSTALL_INTERFACE:include>
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/rulemonkey
)

# Stable namespaced alias for in-tree consumers (`add_subdirectory`):
#   target_link_libraries(my_app PRIVATE RuleMonkey::rulemonkey)
# Mirrors the alias produced by find_package() against an install.
add_library(RuleMonkey::rulemonkey ALIAS rulemonkey)

target_compile_features(rulemonkey PUBLIC cxx_std_17)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(rulemonkey PRIVATE
        -Wall -Wextra -Wno-deprecated-declarations
    )
endif()

# ---------------------------------------------------------------------------
# Optional dev-time profiler instrumentation.
#
# When ON, defines RM_DEV_PROFILES so engine_profile.hpp's gates resolve to
# `true`, and the engine emits per-phase profile reports to stderr at the
# end of each run.  Default OFF: profiler scaffolding compiles to nothing
# (every `if constexpr (k*Profile)` block is dead-stripped).
# ---------------------------------------------------------------------------
option(RULEMONKEY_ENABLE_DEV_PROFILES "Compile in engine dev profiler instrumentation" OFF)
if(RULEMONKEY_ENABLE_DEV_PROFILES)
    target_compile_definitions(rulemonkey PRIVATE RM_DEV_PROFILES)
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
# Examples (off by default — host projects opt in)
# ---------------------------------------------------------------------------
option(RULEMONKEY_BUILD_EXAMPLES "Build the C++ embedding examples" OFF)
if(RULEMONKEY_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
enable_testing()
add_subdirectory(tests)

# ---------------------------------------------------------------------------
# Install rules
#
# Layout under CMAKE_INSTALL_PREFIX (default /usr/local on Linux/macOS):
#
#   include/rulemonkey/{simulator,types}.hpp
#   lib/librulemonkey.a
#   bin/rm_driver
#   lib/cmake/RuleMonkey/RuleMonkeyConfig.cmake
#   lib/cmake/RuleMonkey/RuleMonkeyConfigVersion.cmake
#   lib/cmake/RuleMonkey/RuleMonkeyTargets.cmake
#
# Host projects (e.g. BNGsim) consume RM via:
#
#   find_package(RuleMonkey 3.1 CONFIG REQUIRED)
#   target_link_libraries(my_app PRIVATE RuleMonkey::rulemonkey)
#
# Set RULEMONKEY_INSTALL=OFF when RM is consumed via add_subdirectory()
# from a parent project that doesn't want to advertise install rules.
# ---------------------------------------------------------------------------
option(RULEMONKEY_INSTALL "Generate install() rules for RuleMonkey" ON)

if(RULEMONKEY_INSTALL)
    include(GNUInstallDirs)
    include(CMakePackageConfigHelpers)

    # Public targets: the static library and the rm_driver CLI.
    install(TARGETS rulemonkey rm_driver
        EXPORT RuleMonkeyTargets
        ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

    # Public headers: simulator.hpp + types.hpp under include/rulemonkey/.
    install(DIRECTORY include/rulemonkey
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.hpp"
    )

    # Generate and install the package config files so downstream
    # `find_package(RuleMonkey CONFIG REQUIRED)` works.
    set(_rulemonkey_cmake_dir ${CMAKE_INSTALL_LIBDIR}/cmake/RuleMonkey)

    install(EXPORT RuleMonkeyTargets
        FILE RuleMonkeyTargets.cmake
        NAMESPACE RuleMonkey::
        DESTINATION ${_rulemonkey_cmake_dir}
    )

    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/RuleMonkeyConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
    )

    configure_package_config_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/RuleMonkeyConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/RuleMonkeyConfig.cmake"
        INSTALL_DESTINATION ${_rulemonkey_cmake_dir}
    )

    install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/RuleMonkeyConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/RuleMonkeyConfigVersion.cmake"
        DESTINATION ${_rulemonkey_cmake_dir}
    )
endif()
