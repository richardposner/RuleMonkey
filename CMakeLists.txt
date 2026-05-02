cmake_minimum_required(VERSION 3.25)

project(RuleMonkey VERSION 3.1.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(RM_UNIVERSAL_BINARY
    "On macOS, build a universal arm64+x86_64 binary (~2x build time)" OFF)
if(APPLE AND RM_UNIVERSAL_BINARY)
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

# -Werror is on by default for the in-tree build so a warning shows up
# as a compile failure on the developer's machine before it lands in
# CI.  Downstream consumers building RM as a subdirectory or against an
# installed package can opt out with -DRULEMONKEY_WARNINGS_AS_ERRORS=OFF
# (their toolchain may flag things ours does not — e.g., a newer GCC's
# -Wmaybe-uninitialized tightening).
option(RULEMONKEY_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" ON)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(rulemonkey PRIVATE
        -Wall -Wextra -Wno-deprecated-declarations
    )
    if(RULEMONKEY_WARNINGS_AS_ERRORS)
        target_compile_options(rulemonkey PRIVATE -Werror)
    endif()
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
    target_compile_options(rm_driver PRIVATE -Wall -Wextra -Wno-deprecated-declarations)
    if(RULEMONKEY_WARNINGS_AS_ERRORS)
        target_compile_options(rm_driver PRIVATE -Werror)
    endif()
endif()

# ---------------------------------------------------------------------------
# Optional AddressSanitizer + UndefinedBehaviorSanitizer.
#
# Heavy index-based pool manipulation, BFS / Fenwick code, and pattern-matching
# fast paths are exactly the surface where ASan + UBSan routinely catch
# regressions that observable-trajectory parity tests miss.  Off by default
# (the perf hit is ~2x); cmake --preset asan turns it on.
#
# Both Linux and macOS are exercised in CI (the asan job's matrix).  The
# sanitizer flags themselves are platform-neutral; the CI step sets
# ASAN_OPTIONS / UBSAN_OPTIONS to keep diagnostic output uniform.
# Applied after rm_driver is declared so both library and CLI get the flags.
# ---------------------------------------------------------------------------
option(RULEMONKEY_ENABLE_ASAN "Build with AddressSanitizer + UBSan" OFF)
if(RULEMONKEY_ENABLE_ASAN)
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        message(FATAL_ERROR "RULEMONKEY_ENABLE_ASAN requires Clang or GCC")
    endif()
    set(_rm_asan_flags -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_compile_options(rulemonkey  PRIVATE ${_rm_asan_flags})
    target_link_options   (rulemonkey  PUBLIC  ${_rm_asan_flags})
    target_compile_options(rm_driver   PRIVATE ${_rm_asan_flags})
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

if(RULEMONKEY_INSTALL AND RULEMONKEY_ENABLE_ASAN)
    # Guard against installing an asan-instrumented archive: the asan flags
    # are PRIVATE compile / PUBLIC link on the rulemonkey target, so every
    # downstream `find_package(RuleMonkey)` consumer would inherit
    # `-fsanitize=address,undefined` on its final link line.  That fails
    # outright unless the consumer also enables asan, and even when it
    # links it produces a binary that runs against the wrong runtime.
    # Asan builds are dev-only; the supported `asan` preset disables
    # install for this reason.  Refuse loudly rather than silently
    # producing a poisoned install tree.
    message(FATAL_ERROR
        "RULEMONKEY_ENABLE_ASAN=ON is incompatible with RULEMONKEY_INSTALL=ON: "
        "the asan runtime would propagate to every consumer via PUBLIC link "
        "options and break their builds.  Pass -DRULEMONKEY_INSTALL=OFF (or "
        "use the `asan` preset, which does this for you).")
endif()

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
