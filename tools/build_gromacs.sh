#!/usr/bin/env bash
# SDD:specs.md:§3.4 — Build GROMACS 2027.0-dev from source.
#
# Usage: bash tools/build_gromacs.sh [REPO_ROOT]
#
# Builds GROMACS out-of-tree in external/gromacs/_build and installs to
# external/gromacs/_install. Requires cmake and g++ (C++17).

set -euo pipefail

REPO_ROOT="${1:-$(cd "$(dirname "$0")/.." && pwd)}"
SRC_DIR="$REPO_ROOT/external/gromacs"
BUILD_DIR="$SRC_DIR/_build"
INSTALL_DIR="$SRC_DIR/_install"

# --- Preflight checks ---
CC="${CC:-$(command -v gcc-13 || command -v gcc)}"
CXX="${CXX:-$(command -v g++-13 || command -v g++)}"

for tool in cmake "$CC" "$CXX"; do
    if [[ -z "$tool" ]] || ! command -v "$tool" &>/dev/null; then
        echo "ERROR: Required tool not found: $tool" >&2
        exit 1
    fi
done

echo "=== GROMACS Build Configuration ==="
echo "Source:  $SRC_DIR"
echo "Build:   $BUILD_DIR"
echo "Install: $INSTALL_DIR"
echo "CC:      $CC"
echo "CXX:     $CXX"
echo ""

# --- Configure ---
echo "--- Configuring GROMACS ---"
cmake -B "$BUILD_DIR" -S "$SRC_DIR" \
    -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="$CC" \
    -DCMAKE_CXX_COMPILER="$CXX" \
    -DGMX_MPI=OFF \
    -DGMX_OPENMP=ON \
    -DGMX_GPU=OFF \
    -DGMX_DOUBLE=OFF \
    -DGMX_FFT_LIBRARY=fftpack \
    -DBUILD_TESTING=OFF \
    -DGMX_BUILD_HELP=OFF \
    -DGMX_INSTALL_NBLIB_API=OFF
echo ""

# --- Build ---
echo "--- Building GROMACS ($(nproc) parallel jobs) ---"
cmake --build "$BUILD_DIR" -j"$(nproc)"
echo ""

# --- Install ---
echo "--- Installing GROMACS ---"
cmake --install "$BUILD_DIR"
echo ""

echo "=== GROMACS build complete ==="
echo "Install prefix: $INSTALL_DIR"
