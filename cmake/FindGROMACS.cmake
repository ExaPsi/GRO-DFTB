# SDD:specs.md:§3.2 — Find module for GROMACS shared library and headers.
#
# Searches for libgromacs using GROMACS's own CMake config first (config mode),
# then falls back to manual library/header search.
# Honors GROMACS_ROOT as a search hint (specs.md §3.2).
#
# Defines:
#   GROMACS_FOUND          - True if GROMACS library found
#   GROMACS_INCLUDE_DIRS   - Include directories
#   GROMACS_LIBRARIES      - Library path(s)
#
# Imported target:
#   GROMACS::GROMACS       - Shared library with include dirs set
#
# Note: This find module is for M0/M1. Full GROMACS integration
# (GMX_DFTB flag, DFTBForceProvider wiring) is implemented in US-025 (M2a).

include(FindPackageHandleStandardArgs)

# --- Strategy 1: Use GROMACS's own CMake config ---
# GROMACS installs gromacs-config.cmake under share/cmake/gromacs/.
find_package(GROMACS QUIET CONFIG
    HINTS
        ${GROMACS_ROOT}
        $ENV{GROMACS_ROOT}
    PATH_SUFFIXES
        share/cmake/gromacs
        lib/cmake/gromacs
)

if(GROMACS_FOUND)
    # GROMACS config mode sets Gromacs::libgromacs target.
    # Re-export as GROMACS::GROMACS for our project's conventions.
    if(TARGET Gromacs::libgromacs AND NOT TARGET GROMACS::GROMACS)
        add_library(GROMACS::GROMACS ALIAS Gromacs::libgromacs)
    endif()
    # Extract library and include info from target if available
    if(TARGET Gromacs::libgromacs)
        get_target_property(_gmx_loc Gromacs::libgromacs IMPORTED_LOCATION_RELEASE)
        if(NOT _gmx_loc)
            get_target_property(_gmx_loc Gromacs::libgromacs IMPORTED_LOCATION)
        endif()
        if(_gmx_loc)
            set(GROMACS_LIBRARIES "${_gmx_loc}")
        endif()
        get_target_property(_gmx_inc Gromacs::libgromacs INTERFACE_INCLUDE_DIRECTORIES)
        if(_gmx_inc)
            set(GROMACS_INCLUDE_DIRS "${_gmx_inc}")
        endif()
    endif()
    message(STATUS "Found GROMACS (config mode): ${GROMACS_LIBRARIES}")
    return()
endif()

# --- Strategy 2: Manual search (fallback) ---
find_library(GROMACS_LIBRARY
    NAMES gromacs gromacs_mpi
    HINTS
        ${GROMACS_ROOT}
        $ENV{GROMACS_ROOT}
    PATH_SUFFIXES lib lib64
)

# GROMACS may install headers under gromacs/ or gmxapi/; check both.
find_path(GROMACS_INCLUDE_DIR
    NAMES gromacs/version.h gmxapi/version.h
    HINTS
        ${GROMACS_ROOT}
        $ENV{GROMACS_ROOT}
    PATH_SUFFIXES include
)

find_package_handle_standard_args(GROMACS
    REQUIRED_VARS GROMACS_LIBRARY GROMACS_INCLUDE_DIR
)

if(GROMACS_FOUND)
    set(GROMACS_LIBRARIES ${GROMACS_LIBRARY})
    set(GROMACS_INCLUDE_DIRS ${GROMACS_INCLUDE_DIR})

    if(NOT TARGET GROMACS::GROMACS)
        add_library(GROMACS::GROMACS SHARED IMPORTED)
        set_target_properties(GROMACS::GROMACS PROPERTIES
            IMPORTED_LOCATION "${GROMACS_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${GROMACS_INCLUDE_DIR}"
        )
    endif()
    message(STATUS "Found GROMACS (manual): ${GROMACS_LIBRARIES}")
endif()

mark_as_advanced(GROMACS_LIBRARY GROMACS_INCLUDE_DIR)
