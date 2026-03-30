#!/usr/bin/env bash
# SDD:specs.md:§3 — Verify CMake build system configuration.
#
# Usage: bash tests/test_cmake.sh [REPO_ROOT]
# Exit code: 0 if all tests pass, 1 if any fail.

set -uo pipefail

REPO_ROOT="${1:-$(cd "$(dirname "$0")/.." && pwd)}"
BUILD_DIR=$(mktemp -d /tmp/grodftb_cmake_test.XXXXXX)
LOG_FILE="$BUILD_DIR/cmake_output.log"
PASS=0
FAIL=0

pass() { PASS=$((PASS + 1)); echo "  PASS: $1"; }
fail() { FAIL=$((FAIL + 1)); echo "  FAIL: $1"; }

cleanup() { rm -rf "$BUILD_DIR"; }
trap cleanup EXIT

# ---------------------------------------------------------------------------
# test_cmake_configures — AC-1: cmake configures without errors
# ---------------------------------------------------------------------------
test_cmake_configures() {
    echo "=== test_cmake_configures ==="
    rm -rf "$BUILD_DIR"
    mkdir -p "$BUILD_DIR"
    if cmake -B "$BUILD_DIR" -S "$REPO_ROOT" \
        -DGRODFTB_WITH_GROMACS=OFF \
        -DCMAKE_BUILD_TYPE=Release \
        > "$LOG_FILE" 2>&1; then
        pass "cmake configures successfully"
    else
        fail "cmake configuration failed"
        tail -20 "$LOG_FILE" 2>/dev/null
    fi
}

# ---------------------------------------------------------------------------
# test_cmake_options — AC-3: CMake options match specs.md §3.5
# ---------------------------------------------------------------------------
test_cmake_options() {
    echo "=== test_cmake_options ==="
    if [[ ! -f "$BUILD_DIR/CMakeCache.txt" ]]; then
        fail "CMakeCache.txt not found — run test_cmake_configures first"
        return
    fi

    local cache="$BUILD_DIR/CMakeCache.txt"

    # Check each option and its default
    local -A expected=(
        ["GRODFTB_BUILD_TESTS:BOOL"]="ON"
        ["GRODFTB_BUILD_BENCHMARKS:BOOL"]="OFF"
        ["GRODFTB_WITH_GROMACS:BOOL"]="OFF"
        ["GRODFTB_WITH_PLUMED:BOOL"]="OFF"
        ["GRODFTB_WITH_INAQS:BOOL"]="OFF"
        ["GRODFTB_DOUBLE_PRECISION:BOOL"]="ON"
    )

    for key in "${!expected[@]}"; do
        local var_name="${key%%:*}"
        local val="${expected[$key]}"
        if grep -q "^${key}=${val}$" "$cache"; then
            pass "option $var_name = $val"
        elif grep -q "^${var_name}" "$cache"; then
            local actual
            actual=$(grep "^${key}=" "$cache" | cut -d= -f2)
            # GRODFTB_WITH_GROMACS was explicitly set to OFF in our configure
            if [[ "$var_name" == "GRODFTB_WITH_GROMACS" && "$actual" == "OFF" ]]; then
                pass "option $var_name = $actual (explicitly set)"
            else
                fail "option $var_name expected $val, got $actual"
            fi
        else
            fail "option $var_name not found in cache"
        fi
    done
}

# ---------------------------------------------------------------------------
# test_cmake_version — AC-4: CMakeLists.txt requires minimum version
# ---------------------------------------------------------------------------
test_cmake_version() {
    echo "=== test_cmake_version ==="
    if grep -q "cmake_minimum_required" "$REPO_ROOT/CMakeLists.txt"; then
        pass "cmake_minimum_required declared"
    else
        fail "cmake_minimum_required not found"
    fi
}

# ---------------------------------------------------------------------------
# test_cmake_standards — AC-5: C11 and C++17 standards
# ---------------------------------------------------------------------------
test_cmake_standards() {
    echo "=== test_cmake_standards ==="
    if [[ ! -f "$BUILD_DIR/CMakeCache.txt" ]]; then
        fail "CMakeCache.txt not found"
        return
    fi

    local cache="$BUILD_DIR/CMakeCache.txt"

    if grep -q "CMAKE_C_STANDARD:STRING=11" "$cache"; then
        pass "C_STANDARD = 11"
    else
        fail "C_STANDARD != 11"
    fi

    if grep -q "CMAKE_CXX_STANDARD:STRING=17" "$cache"; then
        pass "CXX_STANDARD = 17"
    else
        fail "CXX_STANDARD != 17"
    fi
}

# ---------------------------------------------------------------------------
# test_find_dftbplus — AC-2: FindDFTBPlus.cmake exists and is used
# ---------------------------------------------------------------------------
test_find_dftbplus() {
    echo "=== test_find_dftbplus ==="
    if [[ -f "$REPO_ROOT/cmake/FindDFTBPlus.cmake" ]]; then
        pass "FindDFTBPlus.cmake exists"
    else
        fail "FindDFTBPlus.cmake missing"
        return
    fi

    # Without a DFTB+ install, verify it handles not-found gracefully
    if [[ -f "$LOG_FILE" ]] && ! grep -qi "fatal\|error.*find.*dftb" "$LOG_FILE"; then
        pass "FindDFTBPlus handles not-found gracefully"
    else
        fail "FindDFTBPlus error on not-found"
    fi
}

# ---------------------------------------------------------------------------
# test_cmake_target_exists — AC-6: libgrodftb target defined
# ---------------------------------------------------------------------------
test_cmake_target_exists() {
    echo "=== test_cmake_target_exists ==="
    if [[ ! -f "$BUILD_DIR/CMakeCache.txt" ]]; then
        fail "CMakeCache.txt not found"
        return
    fi

    # Try to build the grodftb target
    if cmake --build "$BUILD_DIR" --target grodftb > /dev/null 2>&1; then
        pass "grodftb target builds"
    else
        fail "grodftb target does not build"
    fi
}

# ---------------------------------------------------------------------------
# Run all tests
# ---------------------------------------------------------------------------
echo "CMake build system tests for: $REPO_ROOT"
echo ""

test_cmake_configures
echo ""
test_cmake_options
echo ""
test_cmake_version
echo ""
test_cmake_standards
echo ""
test_find_dftbplus
echo ""
test_cmake_target_exists
echo ""

echo "=== Summary: $PASS passed, $FAIL failed ==="
if [[ $FAIL -gt 0 ]]; then
    exit 1
fi
exit 0
