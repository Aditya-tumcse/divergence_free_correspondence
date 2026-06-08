cmake_minimum_required(VERSION 3.20)

project(
    divergence_free_correspondence
    VERSION 0.1.0
    LANGUAGES CXX
)

# ---------------------------------------------------------------------------
# Guard against in-source builds
# ---------------------------------------------------------------------------
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR
        "In-source builds are not allowed.\n"
        "Run: cmake -S . -B build --preset release"
    )
endif()

# ---------------------------------------------------------------------------
# Global C++ settings
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release
        CACHE STRING "Build type (Debug|Release|RelWithDebInfo|MinSizeRel)" FORCE)
endif()

# ---------------------------------------------------------------------------
# Compiler options INTERFACE target — linked by all targets in this project
# ---------------------------------------------------------------------------
add_library(project_options INTERFACE)

target_compile_options(project_options INTERFACE
    $<$<CXX_COMPILER_ID:GNU,Clang>:
        -Wall
        -Wextra
        -Wpedantic
        -Wshadow
        -Wconversion
        -Wsign-conversion
        -Wdouble-promotion
        -Wnon-virtual-dtor
        -Wold-style-cast
        -Woverloaded-virtual
        -Wunused
        -Wno-unused-parameter
    >
    $<$<AND:$<CXX_COMPILER_ID:GNU,Clang>,$<CONFIG:Release>>:
        -O3
        -march=native
        -DNDEBUG
    >
    $<$<AND:$<CXX_COMPILER_ID:GNU,Clang>,$<CONFIG:RelWithDebInfo>>:
        -O2
        -g
        -DNDEBUG
    >
    $<$<AND:$<CXX_COMPILER_ID:GNU,Clang>,$<CONFIG:Debug>>:
        -O0
        -g3
        -fno-omit-frame-pointer
    >
)

# ---------------------------------------------------------------------------
# Dependencies (resolved by Conan-generated toolchain)
# ---------------------------------------------------------------------------
find_package(PCL REQUIRED)
find_package(Ceres REQUIRED)

# ---------------------------------------------------------------------------
# Sub-projects
# ---------------------------------------------------------------------------
add_subdirectory(divergence_free_corr)

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
enable_testing()
add_subdirectory(tests)

# ---------------------------------------------------------------------------
# clang-format targets
#   cmake --build build --target format        — apply formatting in-place
#   cmake --build build --target format-check  — CI-safe dry-run (fails on diff)
# ---------------------------------------------------------------------------
find_program(CLANG_FORMAT clang-format)
if(CLANG_FORMAT)
    file(GLOB_RECURSE ALL_CXX_SOURCES
        "${CMAKE_SOURCE_DIR}/divergence_free_corr/include/*.hpp"
        "${CMAKE_SOURCE_DIR}/divergence_free_corr/src/*.cpp"
        "${CMAKE_SOURCE_DIR}/divergence_free_corr/src/*.hpp"
    )

    add_custom_target(format
        COMMAND ${CLANG_FORMAT} -i ${ALL_CXX_SOURCES}
        COMMENT "clang-format: applying formatting in-place"
        VERBATIM
    )

    add_custom_target(format-check
        COMMAND ${CLANG_FORMAT} --dry-run --Werror ${ALL_CXX_SOURCES}
        COMMENT "clang-format: checking for violations (dry-run)"
        VERBATIM
    )
else()
    message(WARNING "clang-format not found — format and format-check targets unavailable")
endif()

# ---------------------------------------------------------------------------
# cppcheck target
#   cmake --build build --target cppcheck
# ---------------------------------------------------------------------------
find_program(CPPCHECK cppcheck)
if(CPPCHECK)
    add_custom_target(cppcheck
        COMMAND ${CPPCHECK}
            --enable=all
            --std=c++20
            --suppress=missingInclude
            --suppress=missingIncludeSystem
            --suppress=unmatchedSuppression
            --error-exitcode=1
            ${CMAKE_SOURCE_DIR}/divergence_free_corr/src/
        COMMENT "cppcheck: static analysis"
        VERBATIM
    )
else()
    message(STATUS "cppcheck not found — cppcheck target unavailable")
endif()
