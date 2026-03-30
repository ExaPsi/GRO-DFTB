# SDD:specs.md:§3.2 — Find module for DFTB+ shared library and C API header.
#
# Searches for libdftbplus and the dftbplus.h C API header.
# Honors DFTBPLUS_ROOT as a search hint (specs.md §3.2).
#
# Defines:
#   DFTBPlus_FOUND          - True if DFTB+ library and header found
#   DFTBPlus_INCLUDE_DIRS   - Include directories for dftbplus.h
#   DFTBPlus_LIBRARIES      - Library path(s)
#
# Imported target:
#   DFTBPlus::DFTBPlus      - Shared library with include dirs set

include(FindPackageHandleStandardArgs)

# Search for the shared library
find_library(DFTBPlus_LIBRARY
    NAMES dftbplus
    HINTS
        ${DFTBPLUS_ROOT}
        $ENV{DFTBPLUS_ROOT}
    PATH_SUFFIXES lib lib64
)

# Search for the C API header
find_path(DFTBPlus_INCLUDE_DIR
    NAMES dftbplus.h
    HINTS
        ${DFTBPLUS_ROOT}
        $ENV{DFTBPLUS_ROOT}
    PATH_SUFFIXES include include/dftbplus
)

find_package_handle_standard_args(DFTBPlus
    REQUIRED_VARS DFTBPlus_LIBRARY DFTBPlus_INCLUDE_DIR
)

if(DFTBPlus_FOUND)
    set(DFTBPlus_LIBRARIES ${DFTBPlus_LIBRARY})
    set(DFTBPlus_INCLUDE_DIRS ${DFTBPlus_INCLUDE_DIR})

    if(NOT TARGET DFTBPlus::DFTBPlus)
        add_library(DFTBPlus::DFTBPlus SHARED IMPORTED)
        set_target_properties(DFTBPlus::DFTBPlus PROPERTIES
            IMPORTED_LOCATION "${DFTBPlus_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${DFTBPlus_INCLUDE_DIR}"
        )
    endif()
endif()

mark_as_advanced(DFTBPlus_LIBRARY DFTBPlus_INCLUDE_DIR)
