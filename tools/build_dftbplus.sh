#!/usr/bin/env bash
# SDD:specs.md:§3.4 — Build DFTB+ 25.1 from source with API and shared libs.
#
# Usage: bash tools/build_dftbplus.sh [REPO_ROOT]
#
# Builds DFTB+ out-of-tree in external/dftbplus/_build and installs to
# external/dftbplus/_install. Requires cmake, gcc, and gfortran.

set -euo pipefail

REPO_ROOT="${1:-$(cd "$(dirname "$0")/.." && pwd)}"
SRC_DIR="$REPO_ROOT/external/dftbplus"
BUILD_DIR="$SRC_DIR/_build"
INSTALL_DIR="$SRC_DIR/_install"

# --- Preflight checks ---
FC="${FC:-$(command -v gfortran-13 || command -v gfortran)}"
CC="${CC:-$(command -v gcc-13 || command -v gcc)}"
CXX="${CXX:-$(command -v g++-13 || command -v g++)}"

for tool in cmake "$FC" "$CC"; do
    if [[ -z "$tool" ]] || ! command -v "$tool" &>/dev/null; then
        echo "ERROR: Required tool not found: $tool" >&2
        exit 1
    fi
done

echo "=== DFTB+ Build Configuration ==="
echo "Source:  $SRC_DIR"
echo "Build:   $BUILD_DIR"
echo "Install: $INSTALL_DIR"
echo ""

# --- Initialize submodules if needed ---
if [[ -f "$SRC_DIR/.gitmodules" ]]; then
    echo "--- Initializing DFTB+ submodules ---"
    (cd "$SRC_DIR" && git submodule update --init --recursive)
    echo ""
fi

# --- Configure ---
echo "--- Configuring DFTB+ ---"
cmake -B "$BUILD_DIR" -S "$SRC_DIR" \
    -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_Fortran_COMPILER="$FC" \
    -DCMAKE_C_COMPILER="$CC" \
    -DCMAKE_CXX_COMPILER="$CXX" \
    -DWITH_API=ON \
    -DBUILD_SHARED_LIBS=ON \
    -DWITH_OMP=OFF \
    -DWITH_MPI=OFF \
    -DWITH_TBLITE=OFF \
    -DWITH_SDFTD3=OFF
echo ""

# --- Build ---
echo "--- Building DFTB+ ($(nproc) parallel jobs) ---"
cmake --build "$BUILD_DIR" -j"$(nproc)"
echo ""

# --- Install ---
echo "--- Installing DFTB+ ---"
cmake --install "$BUILD_DIR"
echo ""

echo "=== DFTB+ build complete ==="
echo "Install prefix: $INSTALL_DIR"
