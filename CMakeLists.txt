# SDD:specs.md:§3.1 — Top-level build configuration for GRO-DFTB
#
# GRO-DFTB: GROMACS-DFTB+ QM/MM coupling layer
# License: LGPL-3.0-or-later

cmake_minimum_required(VERSION 3.28)

project(gro-dftb
    VERSION 0.2.0
    DESCRIPTION "QM/MM coupling layer integrating DFTB+ with GROMACS"
    LANGUAGES C CXX
)

# US-024: Extract git commit hash for grodftb_version()
execute_process(
    COMMAND git rev-parse --short HEAD
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    OUTPUT_VARIABLE GRODFTB_GIT_HASH
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT GRODFTB_GIT_HASH)
    set(GRODFTB_GIT_HASH "unknown")
endif()
configure_file(
    "${CMAKE_SOURCE_DIR}/src/driver/version.c.in"
    "${CMAKE_CURRENT_BINARY_DIR}/version.c"
    @ONLY
)

# SDD:specs.md:§3.3 — Compiler standards
set(CMAKE_C_STANDARD 11 CACHE STRING "C standard version")
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard version")
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# SDD:specs.md:§3.5 — Project options
option(GRODFTB_BUILD_TESTS      "Build test suite"              ON)
option(GRODFTB_BUILD_BENCHMARKS "Build benchmark runner"        OFF)
option(GRODFTB_WITH_GROMACS     "Build GROMACS backend"         ON)
set(GRODFTB_GROMACS_ROUTE "A" CACHE STRING "GROMACS integration route (A or B)")
option(GRODFTB_WITH_PLUMED      "Build PLUMED CV plugin"        OFF)
option(GRODFTB_WITH_INAQS       "Build INAQS backend"           OFF)
option(GRODFTB_DOUBLE_PRECISION "Use double precision throughout" ON)

# US-037: Test gating options
option(GRODFTB_LONG_TESTS       "Build long-running tests (100+ ps)" OFF)
option(GRODFTB_MEMORY_TESTS     "Build memory safety tests (Valgrind/ASan)" OFF)

# Fortran is needed for linking DFTB+; enable it only if a compiler is found.
include(CheckLanguage)
check_language(Fortran)
if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
else()
    message(STATUS "No Fortran compiler found — DFTB+ linking will require one at M1")
endif()

# Module path for custom find modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# ---------------------------------------------------------------------------
# Find DFTB+ (optional at M0; required from M1 onward)
# ---------------------------------------------------------------------------
# Default location: external/dftbplus/_install (per CLAUDE.md)
if(NOT DFTBPLUS_ROOT)
    set(_dftbplus_default "${CMAKE_SOURCE_DIR}/external/dftbplus/_install")
    if(EXISTS "${_dftbplus_default}/lib/libdftbplus.so")
        set(DFTBPLUS_ROOT "${_dftbplus_default}" CACHE PATH "Path to DFTB+ install")
    endif()
endif()

find_package(DFTBPlus QUIET)
if(DFTBPlus_FOUND)
    message(STATUS "Found DFTB+: ${DFTBPlus_LIBRARIES}")
else()
    message(STATUS "DFTB+ not found — libgrodftb will be built without DFTB+ linkage")
endif()

# ---------------------------------------------------------------------------
# Find GROMACS (optional at M0; required from M2 onward)
# ---------------------------------------------------------------------------
if(GRODFTB_WITH_GROMACS)
    # Default location: external/gromacs/_install (per CLAUDE.md)
    if(NOT GROMACS_DIR)
        set(_gromacs_default "${CMAKE_SOURCE_DIR}/external/gromacs/_install/share/cmake/gromacs")
        if(EXISTS "${_gromacs_default}/gromacs-config.cmake")
            set(GROMACS_DIR "${_gromacs_default}" CACHE PATH "Path to GROMACS CMake config")
        endif()
    endif()

    find_package(GROMACS QUIET)
    if(GROMACS_FOUND)
        message(STATUS "Found GROMACS: ${GROMACS_LIBRARIES}")
    else()
        message(STATUS "GROMACS not found — GROMACS backend will not be built")
        message(STATUS "  Hint: Set GROMACS_DIR to the GROMACS cmake config directory")
    endif()
endif()

# ---------------------------------------------------------------------------
# libgrodftb shared library (SDD:specs.md:§3.1)
# ---------------------------------------------------------------------------
add_library(grodftb SHARED
    src/driver/dftb_driver.c
    src/units/units.c
    src/linkatom/linkatom.c
    src/embedding/switching.c
    src/embedding/damping.c
    ${CMAKE_CURRENT_BINARY_DIR}/version.c
)

target_include_directories(grodftb PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(grodftb PRIVATE m)

if(DFTBPlus_FOUND)
    target_link_libraries(grodftb PRIVATE DFTBPlus::DFTBPlus)
    target_compile_definitions(grodftb PRIVATE GRODFTB_HAS_DFTBPLUS=1)
endif()

set_target_properties(grodftb PROPERTIES
    OUTPUT_NAME grodftb
    VERSION ${PROJECT_VERSION}
    SOVERSION 0
)

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
if(GRODFTB_BUILD_TESTS)
    enable_testing()
    # US-006: Unit conversion tests
    add_executable(test_units tests/test_units.c)
    target_link_libraries(test_units PRIVATE grodftb m)
    add_test(NAME test_units COMMAND test_units)

    # US-009: Driver init/finalize tests
    add_executable(test_driver tests/test_driver.c)
    target_link_libraries(test_driver PRIVATE grodftb m)
    target_compile_definitions(test_driver PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_driver PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_driver PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    add_test(NAME test_driver COMMAND test_driver)

    # US-015: Finite-difference force validation
    add_executable(test_forces tests/test_forces.c)
    target_link_libraries(test_forces PRIVATE grodftb m)
    target_compile_definitions(test_forces PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_forces PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_forces PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    add_test(NAME test_forces COMMAND test_forces)
    set_tests_properties(test_forces PROPERTIES LABELS "finitediff")

    # US-033: Link atom placement tests
    add_executable(test_linkatom tests/test_linkatom.c)
    target_link_libraries(test_linkatom PRIVATE grodftb m)
    add_test(NAME test_linkatom COMMAND test_linkatom)

    # US-034: Link atom force projection FD validation
    add_executable(test_linkatom_fd tests/test_linkatom_fd.c)
    target_link_libraries(test_linkatom_fd PRIVATE grodftb m)
    target_compile_definitions(test_linkatom_fd PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_linkatom_fd PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_linkatom_fd PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    add_test(NAME test_linkatom_fd COMMAND test_linkatom_fd)
    set_tests_properties(test_linkatom_fd PROPERTIES LABELS "finitediff")

    # US-035: Charge redistribution for link atom boundaries
    add_executable(test_charge_redistrib tests/test_charge_redistrib.c)
    target_link_libraries(test_charge_redistrib PRIVATE grodftb m)
    target_compile_definitions(test_charge_redistrib PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_charge_redistrib PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_charge_redistrib PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    add_test(NAME test_charge_redistrib COMMAND test_charge_redistrib)
    set_tests_properties(test_charge_redistrib PROPERTIES LABELS "finitediff")

    # US-036: FD force validation for B6 alanine dipeptide
    # ThDD:T-US-036-V.4 — Central difference FD validation
    # SDD:specs.md:§21.2 — B6 benchmark system
    add_executable(test_linkatom_fd_b6 tests/test_linkatom_fd_b6.c)
    target_link_libraries(test_linkatom_fd_b6 PRIVATE grodftb m)
    target_compile_definitions(test_linkatom_fd_b6 PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_linkatom_fd_b6 PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_linkatom_fd_b6 PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    add_test(NAME test_linkatom_fd_b6 COMMAND test_linkatom_fd_b6)
    set_tests_properties(test_linkatom_fd_b6 PROPERTIES LABELS "finitediff;b6")

    # -------------------------------------------------------------------------
    # US-037: Integration tests on B6 (alanine dipeptide)
    # ThDD:T-US-037-V.1 — Force continuity across QM/MM boundary
    # ThDD:T-US-037-N.7 — NVE drift measurement
    # SDD:specs.md:§21.2 — B6 benchmark system
    # -------------------------------------------------------------------------

    # Main integration tests (force continuity, bond stability)
    add_executable(test_b6_integration tests/test_b6_integration.c)
    target_link_libraries(test_b6_integration PRIVATE grodftb m)
    target_compile_definitions(test_b6_integration PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_b6_integration PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_b6_integration PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    add_test(NAME test_b6_integration COMMAND test_b6_integration)
    set_tests_properties(test_b6_integration PROPERTIES LABELS "integration;b6")

    # NVE energy conservation test (long-running, gated)
    if(GRODFTB_LONG_TESTS)
        add_executable(test_b6_nve tests/test_b6_nve.c)
        target_link_libraries(test_b6_nve PRIVATE grodftb m)
        target_compile_definitions(test_b6_nve PRIVATE
            GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
        )
        if(DFTBPlus_FOUND)
            target_link_libraries(test_b6_nve PRIVATE DFTBPlus::DFTBPlus)
            target_compile_definitions(test_b6_nve PRIVATE GRODFTB_HAS_DFTBPLUS=1)
        endif()
        add_test(NAME test_b6_nve COMMAND test_b6_nve)
        set_tests_properties(test_b6_nve PROPERTIES LABELS "nve;b6;long")
    endif()

    # Memory safety wrapper (Valgrind/ASan, gated)
    if(GRODFTB_MEMORY_TESTS)
        add_executable(test_b6_memory tests/test_b6_memory.c)
        target_link_libraries(test_b6_memory PRIVATE grodftb m)
        target_compile_definitions(test_b6_memory PRIVATE
            GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
        )
        if(DFTBPlus_FOUND)
            target_link_libraries(test_b6_memory PRIVATE DFTBPlus::DFTBPlus)
            target_compile_definitions(test_b6_memory PRIVATE GRODFTB_HAS_DFTBPLUS=1)
        endif()
        add_test(NAME test_b6_memory COMMAND test_b6_memory)
        set_tests_properties(test_b6_memory PROPERTIES LABELS "memory;b6")
    endif()

    # -------------------------------------------------------------------------
    # US-038: Cutoff-Based Electrostatic Embedding Tests
    # ThDD:T-US-038-N.1 -- Cutoff neighbor gathering algorithm
    # ThDD:T-US-038-V.1 -- Energy agreement with exact Coulomb sum
    # ThDD:T-US-038-V.2 -- FD force validation
    # SDD:specs.md:S8 -- Module 4: Electrostatic Embedding
    # -------------------------------------------------------------------------

    # NOTE: test_cutoff_embedding.c is DEPRECATED (2026-02-04)
    # It tested a hypothetical standalone libgrodftb API that was never implemented.
    # The actual embedding is in DFTBForceProvider (GROMACS backend).
    # GROMACS-side tests: external/gromacs/src/gromacs/applied_forces/qmmm/tests/
    #                     dftbforceprovider_embedding.cpp
    # The file is retained for documentation but removed from the build.

    # Embedding finite-difference force validation tests (validates actual implementation)
    add_executable(test_embedding_fd tests/test_embedding_fd.c)
    target_link_libraries(test_embedding_fd PRIVATE grodftb m)
    target_compile_definitions(test_embedding_fd PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_embedding_fd PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_embedding_fd PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    add_test(NAME test_embedding_fd COMMAND test_embedding_fd)
    set_tests_properties(test_embedding_fd PROPERTIES LABELS "finitediff;embedding;b4")

    # -------------------------------------------------------------------------
    # US-039: Switching Function for Cutoff Embedding Tests
    # ThDD:T-US-039-V.1 -- Switching function algebraic identities
    # ThDD:T-US-039-V.2 -- First derivative factorization
    # ThDD:T-US-039-V.3 -- Second derivative factorization
    # ThDD:T-US-039-V.4 -- FD force validation in switching region
    # ThDD:T-US-039-V.11 -- Hard cutoff equivalence when switch_width=0
    # ThDD:T-US-039-V.12 -- B4 reference agreement with switching disabled
    # SDD:specs.md:S8.1 -- Switching function specification
    # -------------------------------------------------------------------------

    # Switching function mathematical identity tests (pure function tests)
    add_executable(test_switching tests/test_switching.c)
    target_link_libraries(test_switching PRIVATE grodftb m)
    add_test(NAME test_switching COMMAND test_switching)
    set_tests_properties(test_switching PROPERTIES LABELS "switching;identity")

    # Switching FD force validation tests
    # NOTE: Will fail until grodftb_set_embedding_params() is implemented (TDD RED phase)
    add_executable(test_switching_fd tests/test_switching_fd.c)
    target_link_libraries(test_switching_fd PRIVATE grodftb m)
    target_compile_definitions(test_switching_fd PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_switching_fd PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_switching_fd PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    add_test(NAME test_switching_fd COMMAND test_switching_fd)
    set_tests_properties(test_switching_fd PROPERTIES LABELS "finitediff;switching;b4")

    # Switching regression tests (backward compatibility)
    # NOTE: Will fail until grodftb_set_embedding_params() is implemented (TDD RED phase)
    add_executable(test_switching_regression tests/test_switching_regression.c)
    target_link_libraries(test_switching_regression PRIVATE grodftb m)
    target_compile_definitions(test_switching_regression PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_switching_regression PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_switching_regression PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    add_test(NAME test_switching_regression COMMAND test_switching_regression)
    set_tests_properties(test_switching_regression PROPERTIES LABELS "regression;switching;b4")

    # -------------------------------------------------------------------------
    # US-041: Virial FD Validation Tests
    # ThDD:T-US-041-7.1 -- QM/MM virial tensor from electrostatic embedding
    # ThDD:T-US-041-V.1 -- FD validation via strain perturbation
    # SDD:specs.md:S21.1 -- Finite-difference validation tolerance
    # -------------------------------------------------------------------------
    add_executable(test_virial_fd tests/test_virial_fd.c)
    target_link_libraries(test_virial_fd PRIVATE grodftb m)
    target_compile_definitions(test_virial_fd PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_virial_fd PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_virial_fd PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    add_test(NAME test_virial_fd COMMAND test_virial_fd)
    set_tests_properties(test_virial_fd PROPERTIES LABELS "finitediff;virial;b4")

    # -------------------------------------------------------------------------
    # US-041: B5 Preprocessing Tests
    # ThDD:T-US-041-PRE.1 -- QM topology verification
    # ThDD:T-US-041-PRE.2 -- Species index verification
    # ThDD:T-US-041-PRE.3 -- Embedding charge collection verification
    # SDD:docs/verification/US-041.md:Section 3.2 -- Pre-simulation tests
    # -------------------------------------------------------------------------
    add_executable(test_b5_preprocessing tests/test_b5_preprocessing.c)
    target_link_libraries(test_b5_preprocessing PRIVATE grodftb m)
    target_compile_definitions(test_b5_preprocessing PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_b5_preprocessing PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_b5_preprocessing PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    add_test(NAME test_b5_preprocessing COMMAND test_b5_preprocessing)
    set_tests_properties(test_b5_preprocessing PROPERTIES LABELS "unit;b5;qmmm")

    # -------------------------------------------------------------------------
    # US-041: B5 FD Force Validation Tests
    # ThDD:T-US-041-STEP-2 -- FD validation for QM forces
    # ThDD:T-US-041-STEP-3 -- FD validation for MM forces (inner/switching/outer)
    # SDD:docs/verification/US-041.md:Section 4 -- FD methodology
    # SDD:specs.md:S21.1 -- Finite-difference validation tolerance
    # -------------------------------------------------------------------------
    add_executable(test_b5_fd_forces tests/test_b5_fd_forces.c)
    target_link_libraries(test_b5_fd_forces PRIVATE grodftb m)
    target_compile_definitions(test_b5_fd_forces PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_b5_fd_forces PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_b5_fd_forces PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    add_test(NAME test_b5_fd_forces COMMAND test_b5_fd_forces)
    set_tests_properties(test_b5_fd_forces PROPERTIES LABELS "finitediff;b5;qmmm")

    # -------------------------------------------------------------------------
    # US-041: B5 Single-Step Physics Verification Tests
    # ThDD:T-US-041-2.1 -- Embedding energy consistency
    # ThDD:T-US-041-4.2 -- Newton's third law for QM-MM pairs
    # ThDD:T-US-041-5.3 -- Force continuity at switching boundaries
    # SDD:docs/verification/US-041.md:Section 3.3 -- Single-step tests
    # SDD:specs.md:S21.1 -- Test tolerances
    # -------------------------------------------------------------------------
    add_executable(test_b5_single_step tests/test_b5_single_step.c)
    target_link_libraries(test_b5_single_step PRIVATE grodftb m)
    target_compile_definitions(test_b5_single_step PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_b5_single_step PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_b5_single_step PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    add_test(NAME test_b5_single_step COMMAND test_b5_single_step)
    set_tests_properties(test_b5_single_step PROPERTIES LABELS "integration;b5;qmmm")

    # -------------------------------------------------------------------------
    # US-041: B5 NVE Trajectory Stability Tests
    # ThDD:T-US-041-TRAJ.1 -- NVE energy drift via linear regression
    # ThDD:T-US-041-TRAJ.2 -- Energy fluctuation analysis
    # ThDD:T-US-041-TRAJ.4 -- Numerical error detection
    # ThDD:T-US-041-TRAJ.5 -- Switching region crossing analysis
    # ThDD:T-US-041-TRAJ.6 -- QM charge stability
    # SDD:docs/verification/US-041.md:Section 5 -- Trajectory tests
    # SDD:specs.md:S21.1 -- NVE drift < 0.01 kJ/mol/ps/atom
    # -------------------------------------------------------------------------
    add_executable(test_b5_nve tests/test_b5_nve.c)
    target_link_libraries(test_b5_nve PRIVATE grodftb m)
    target_compile_definitions(test_b5_nve PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_b5_nve PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_b5_nve PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    if(GRODFTB_LONG_TESTS)
        target_compile_definitions(test_b5_nve PRIVATE GRODFTB_LONG_TESTS=1)
    endif()
    add_test(NAME test_b5_nve COMMAND test_b5_nve)
    set_tests_properties(test_b5_nve PROPERTIES LABELS "stability;nve;b5;qmmm")

    # -------------------------------------------------------------------------
    # US-042: NVE Energy Drift Statistical Analysis Tests
    # ThDD:T-US-042-N.1 through N.17 -- Drift measurement methodology
    # ThDD:T-US-042-N.3 -- Linear regression algorithm
    # ThDD:T-US-042-N.4 -- R-squared calculation
    # ThDD:T-US-042-N.6 -- Standard error of slope
    # ThDD:T-US-042-N.9 -- Temperature fluctuation formula
    # SDD:specs.md:S21.1 -- NVE drift criterion (< 0.01 kJ/mol/ps/atom)
    # SDD:docs/verification/US-042.md -- Verification plan
    # -------------------------------------------------------------------------

    # Unit tests: Statistical algorithms (fast, always run)
    add_executable(test_nve_drift tests/test_nve_drift.c)
    target_link_libraries(test_nve_drift PRIVATE grodftb m)
    target_compile_definitions(test_nve_drift PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_nve_drift PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_nve_drift PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    if(GRODFTB_LONG_TESTS)
        target_compile_definitions(test_nve_drift PRIVATE GRODFTB_LONG_TESTS=1)
    endif()

    # Register unit tests (always available, fast)
    add_test(NAME test_nve_drift_unit COMMAND test_nve_drift --unit)
    set_tests_properties(test_nve_drift_unit PROPERTIES
        LABELS "unit;nve;drift;us042"
        TIMEOUT 60
    )

    # Register integration tests (gated, require trajectory data)
    if(GRODFTB_LONG_TESTS)
        add_test(NAME test_nve_drift_b5_10ps COMMAND test_nve_drift --integration)
        set_tests_properties(test_nve_drift_b5_10ps PROPERTIES
            LABELS "integration;nve;drift;b5;us042;long"
            TIMEOUT 14400
        )
    endif()

    # -------------------------------------------------------------------------
    # US-043: NVE Stability Test Suite for B5 QM/MM System
    # ThDD:T-US-043-V.1 through V.4 -- Energy conservation tests
    # ThDD:T-US-043-V.5 through V.7 -- Trajectory validity tests
    # ThDD:T-US-043-V.8 through V.9 -- Temperature tests
    # ThDD:T-US-043-V.10 through V.12 -- QM region stability tests
    # ThDD:T-US-043-V.13 through V.14 -- Checkpoint/restart consistency
    # ThDD:T-US-043-V.15 through V.16 -- Embedding stability tests
    # SDD:specs.md:S21.1 -- NVE drift criterion (< 0.01 kJ/mol/ps/atom)
    # SDD:docs/verification/US-043.md -- Verification plan
    # -------------------------------------------------------------------------

    # Test executable: NVE stability tests (V.1-V.4)
    add_executable(test_nve_stability tests/test_nve_stability.c)
    target_link_libraries(test_nve_stability PRIVATE grodftb m)
    target_compile_definitions(test_nve_stability PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_nve_stability PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_nve_stability PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    if(GRODFTB_LONG_TESTS)
        target_compile_definitions(test_nve_stability PRIVATE GRODFTB_LONG_TESTS=1)
    endif()

    # Test executable: Trajectory validity tests (V.5-V.7)
    add_executable(test_trajectory_validity tests/test_trajectory_validity.c)
    target_link_libraries(test_trajectory_validity PRIVATE grodftb m)
    target_compile_definitions(test_trajectory_validity PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_trajectory_validity PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_trajectory_validity PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    if(GRODFTB_LONG_TESTS)
        target_compile_definitions(test_trajectory_validity PRIVATE GRODFTB_LONG_TESTS=1)
    endif()

    # Test executable: Temperature tests (V.8-V.9)
    add_executable(test_temperature tests/test_temperature.c)
    target_link_libraries(test_temperature PRIVATE grodftb m)
    target_compile_definitions(test_temperature PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_temperature PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_temperature PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    if(GRODFTB_LONG_TESTS)
        target_compile_definitions(test_temperature PRIVATE GRODFTB_LONG_TESTS=1)
    endif()

    # Test executable: QM stability tests (V.10-V.12)
    add_executable(test_qm_stability tests/test_qm_stability.c)
    target_link_libraries(test_qm_stability PRIVATE grodftb m)
    target_compile_definitions(test_qm_stability PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_qm_stability PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_qm_stability PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    if(GRODFTB_LONG_TESTS)
        target_compile_definitions(test_qm_stability PRIVATE GRODFTB_LONG_TESTS=1)
    endif()

    # Test executable: Checkpoint tests (V.13-V.14)
    add_executable(test_checkpoint tests/test_checkpoint.c)
    target_link_libraries(test_checkpoint PRIVATE grodftb m)
    target_compile_definitions(test_checkpoint PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_checkpoint PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_checkpoint PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    if(GRODFTB_LONG_TESTS)
        target_compile_definitions(test_checkpoint PRIVATE GRODFTB_LONG_TESTS=1)
    endif()

    # Test executable: Embedding stability tests (V.15-V.16)
    add_executable(test_embedding_stability tests/test_embedding_stability.c)
    target_link_libraries(test_embedding_stability PRIVATE grodftb m)
    target_compile_definitions(test_embedding_stability PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_embedding_stability PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_embedding_stability PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    if(GRODFTB_LONG_TESTS)
        target_compile_definitions(test_embedding_stability PRIVATE GRODFTB_LONG_TESTS=1)
    endif()

    # Register unit tests (always run, fast)
    add_test(NAME test_nve_stability_unit COMMAND test_nve_stability --unit)
    set_tests_properties(test_nve_stability_unit PROPERTIES
        LABELS "unit;nve;stability;us043"
        TIMEOUT 60
    )

    add_test(NAME test_trajectory_validity_unit COMMAND test_trajectory_validity --unit)
    set_tests_properties(test_trajectory_validity_unit PROPERTIES
        LABELS "unit;nve;stability;us043"
        TIMEOUT 60
    )

    add_test(NAME test_temperature_unit COMMAND test_temperature --unit)
    set_tests_properties(test_temperature_unit PROPERTIES
        LABELS "unit;nve;stability;us043"
        TIMEOUT 60
    )

    add_test(NAME test_qm_stability_unit COMMAND test_qm_stability --unit)
    set_tests_properties(test_qm_stability_unit PROPERTIES
        LABELS "unit;nve;stability;us043"
        TIMEOUT 60
    )

    add_test(NAME test_checkpoint_unit COMMAND test_checkpoint --unit)
    set_tests_properties(test_checkpoint_unit PROPERTIES
        LABELS "unit;nve;stability;us043"
        TIMEOUT 60
    )

    add_test(NAME test_embedding_stability_unit COMMAND test_embedding_stability --unit)
    set_tests_properties(test_embedding_stability_unit PROPERTIES
        LABELS "unit;nve;stability;us043"
        TIMEOUT 60
    )

    # Register integration tests (gated behind GRODFTB_LONG_TESTS, 500ps trajectories)
    if(GRODFTB_LONG_TESTS)
        add_test(NAME test_nve_stability_b5_500ps COMMAND test_nve_stability --integration)
        set_tests_properties(test_nve_stability_b5_500ps PROPERTIES
            LABELS "integration;nve;stability;b5;us043;long"
            TIMEOUT 86400
        )

        add_test(NAME test_trajectory_validity_b5_500ps COMMAND test_trajectory_validity --integration)
        set_tests_properties(test_trajectory_validity_b5_500ps PROPERTIES
            LABELS "integration;nve;stability;b5;us043;long"
            TIMEOUT 86400
        )

        add_test(NAME test_temperature_b5_500ps COMMAND test_temperature --integration)
        set_tests_properties(test_temperature_b5_500ps PROPERTIES
            LABELS "integration;nve;stability;b5;us043;long"
            TIMEOUT 86400
        )

        add_test(NAME test_qm_stability_b5_500ps COMMAND test_qm_stability --integration)
        set_tests_properties(test_qm_stability_b5_500ps PROPERTIES
            LABELS "integration;nve;stability;b5;us043;long"
            TIMEOUT 86400
        )

        add_test(NAME test_checkpoint_b5_500ps COMMAND test_checkpoint --integration)
        set_tests_properties(test_checkpoint_b5_500ps PROPERTIES
            LABELS "integration;nve;stability;b5;us043;long"
            TIMEOUT 86400
        )

        add_test(NAME test_embedding_stability_b5_500ps COMMAND test_embedding_stability --integration)
        set_tests_properties(test_embedding_stability_b5_500ps PROPERTIES
            LABELS "integration;nve;stability;b5;us043;long"
            TIMEOUT 86400
        )
    endif()

    # -------------------------------------------------------------------------
    # US-043b: Gaussian-Damped Embedding Tests
    # ThDD:T-US-043b-1.11 -- Damped Coulomb potential erf(r/sigma)/r
    # ThDD:T-US-043b-2.4  -- Radial derivative
    # ThDD:T-US-043b-9.1  -- Taylor expansion for small r
    # ThDD:T-US-043b-3.6  -- Combined damped+switched kernel
    # ThDD:T-US-043b-4.6  -- MM embedding force
    # SDD:specs.md:S8.1   -- Embedding potential specification
    # SDD:specs.md:S21.1  -- NVE drift criterion
    # -------------------------------------------------------------------------

    # Unit tests: damped Coulomb functions (AC-1..5, AC-15)
    add_executable(test_damping tests/test_damping.c)
    target_link_libraries(test_damping PRIVATE grodftb m)
    target_compile_definitions(test_damping PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_damping PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_damping PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    add_test(NAME test_damping COMMAND test_damping)
    set_tests_properties(test_damping PROPERTIES LABELS "unit;damping;us043b")

    # FD validation tests: damped gradients (AC-6..9)
    add_executable(test_damping_fd tests/test_damping_fd.c)
    target_link_libraries(test_damping_fd PRIVATE grodftb m)
    target_compile_definitions(test_damping_fd PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_damping_fd PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_damping_fd PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    add_test(NAME test_damping_fd COMMAND test_damping_fd)
    set_tests_properties(test_damping_fd PROPERTIES LABELS "finitediff;damping;us043b")

    # Integration tests: B4 damped energy match, QM-QM exclusion (AC-10, AC-14)
    add_executable(test_damping_integration tests/test_damping_integration.c)
    target_link_libraries(test_damping_integration PRIVATE grodftb m)
    target_compile_definitions(test_damping_integration PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_damping_integration PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_damping_integration PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    add_test(NAME test_damping_integration COMMAND test_damping_integration)
    set_tests_properties(test_damping_integration PROPERTIES LABELS "integration;damping;b4;us043b")

    # NVE stability test with damping (AC-11, gated)
    if(GRODFTB_LONG_TESTS)
        add_executable(test_nve_drift_damped tests/test_nve_drift_damped.c)
        target_link_libraries(test_nve_drift_damped PRIVATE grodftb m)
        target_compile_definitions(test_nve_drift_damped PRIVATE
            GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
            GRODFTB_LONG_TESTS=1
        )
        if(DFTBPlus_FOUND)
            target_link_libraries(test_nve_drift_damped PRIVATE DFTBPlus::DFTBPlus)
            target_compile_definitions(test_nve_drift_damped PRIVATE GRODFTB_HAS_DFTBPLUS=1)
        endif()
        add_test(NAME test_nve_drift_damped COMMAND test_nve_drift_damped)
        set_tests_properties(test_nve_drift_damped PROPERTIES
            LABELS "stability;nve;damping;b5;us043b;long"
            TIMEOUT 14400
        )
    endif()

    # ===== US-046: PME potential unit conversion tests (AC-4) =====
    add_executable(test_pme_potential_units tests/test_pme_potential_units.c)
    target_link_libraries(test_pme_potential_units PRIVATE m)
    add_test(NAME test_pme_potential_units COMMAND test_pme_potential_units)
    set_tests_properties(test_pme_potential_units PROPERTIES LABELS "unit;pme;us046")

    # US-046: erfc kernel + damped Ewald tests (AC-7, AC-8, AC-9)
    add_executable(test_pme_erfc_kernel tests/test_pme_erfc_kernel.c)
    target_link_libraries(test_pme_erfc_kernel PRIVATE grodftb m)
    add_test(NAME test_pme_erfc_kernel COMMAND test_pme_erfc_kernel)
    set_tests_properties(test_pme_erfc_kernel PROPERTIES LABELS "unit;pme;erfc;us046")

    # US-048: Embedding energy consistency check
    add_executable(test_consistency_check tests/test_consistency_check.c)
    target_link_libraries(test_consistency_check PRIVATE grodftb m)
    target_compile_definitions(test_consistency_check PRIVATE
        GRODFTB_SOURCE_DIR_STR="${CMAKE_SOURCE_DIR}"
    )
    if(DFTBPlus_FOUND)
        target_link_libraries(test_consistency_check PRIVATE DFTBPlus::DFTBPlus)
        target_compile_definitions(test_consistency_check PRIVATE GRODFTB_HAS_DFTBPLUS=1)
    endif()
    add_test(NAME test_consistency_check COMMAND test_consistency_check)
    set_tests_properties(test_consistency_check PROPERTIES LABELS "integration;consistency;us048")
endif()
