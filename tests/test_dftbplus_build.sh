#!/usr/bin/env bash
# SDD:specs.md:§3.4 — Verify DFTB+ build from source.
#
# Usage: bash tests/test_dftbplus_build.sh [REPO_ROOT]
# Exit code: 0 if all tests pass, 1 if any fail.

set -uo pipefail

REPO_ROOT="${1:-$(cd "$(dirname "$0")/.." && pwd)}"
INSTALL_PREFIX="$REPO_ROOT/external/dftbplus/_install"
PASS=0
FAIL=0

pass() { PASS=$((PASS + 1)); echo "  PASS: $1"; }
fail() { FAIL=$((FAIL + 1)); echo "  FAIL: $1"; }

# ---------------------------------------------------------------------------
# test_build_script_exists — build script exists and is executable
# ---------------------------------------------------------------------------
test_build_script_exists() {
    echo "=== test_build_script_exists ==="
    local script="$REPO_ROOT/tools/build_dftbplus.sh"
    if [[ -f "$script" ]]; then
        pass "tools/build_dftbplus.sh exists"
    else
        fail "tools/build_dftbplus.sh missing"
    fi
}

# ---------------------------------------------------------------------------
# test_dftbplus_lib_exists — libdftbplus.so present under install prefix
# ---------------------------------------------------------------------------
test_dftbplus_lib_exists() {
    echo "=== test_dftbplus_lib_exists ==="
    local found=0
    for d in lib lib64; do
        if [[ -f "$INSTALL_PREFIX/$d/libdftbplus.so" ]]; then
            found=1
            pass "libdftbplus.so found in $INSTALL_PREFIX/$d/"
            break
        fi
    done
    if [[ $found -eq 0 ]]; then
        fail "libdftbplus.so not found under $INSTALL_PREFIX"
    fi
}

# ---------------------------------------------------------------------------
# test_dftbplus_header_exists — dftbplus.h present under install prefix
# ---------------------------------------------------------------------------
test_dftbplus_header_exists() {
    echo "=== test_dftbplus_header_exists ==="
    if find "$INSTALL_PREFIX" -name "dftbplus.h" 2>/dev/null | grep -q .; then
        pass "dftbplus.h found under $INSTALL_PREFIX"
    else
        fail "dftbplus.h not found under $INSTALL_PREFIX"
    fi
}

# ---------------------------------------------------------------------------
# test_dftbplus_api_version — compile and run smoke test, check API 0.4.0
# ---------------------------------------------------------------------------
test_dftbplus_api_version() {
    echo "=== test_dftbplus_api_version ==="
    local smoke_src="$REPO_ROOT/tests/smoke_dftbp_api.c"
    local smoke_bin
    smoke_bin=$(mktemp /tmp/smoke_dftbp_api.XXXXXX)

    if [[ ! -f "$smoke_src" ]]; then
        fail "smoke_dftbp_api.c not found"
        return
    fi

    # Find include and lib dirs
    local inc_dir
    inc_dir=$(find "$INSTALL_PREFIX" -name "dftbplus.h" -printf '%h\n' 2>/dev/null | head -1)
    local lib_dir=""
    for d in lib lib64; do
        if [[ -f "$INSTALL_PREFIX/$d/libdftbplus.so" ]]; then
            lib_dir="$INSTALL_PREFIX/$d"
            break
        fi
    done

    if [[ -z "$inc_dir" || -z "$lib_dir" ]]; then
        fail "cannot locate DFTB+ headers or library for compilation"
        rm -f "$smoke_bin"
        return
    fi

    # Compile
    if ! gcc -o "$smoke_bin" "$smoke_src" -I"$inc_dir" -L"$lib_dir" -ldftbplus -Wl,-rpath,"$lib_dir" 2>/dev/null; then
        fail "smoke_dftbp_api.c failed to compile"
        rm -f "$smoke_bin"
        return
    fi

    # Run
    local output
    output=$("$smoke_bin" 2>&1)
    local rc=$?
    rm -f "$smoke_bin"

    if [[ $rc -ne 0 ]]; then
        fail "smoke test exited with code $rc: $output"
        return
    fi

    if echo "$output" | grep -qE "0[. ]4[. ]0"; then
        pass "API version 0.4.0 confirmed"
    else
        fail "API version mismatch: $output"
    fi
}

# ---------------------------------------------------------------------------
# test_find_dftbplus_with_root — cmake finds DFTB+ via DFTBPLUS_ROOT
# ---------------------------------------------------------------------------
test_find_dftbplus_with_root() {
    echo "=== test_find_dftbplus_with_root ==="
    local build_dir
    build_dir=$(mktemp -d /tmp/grodftb_find_test.XXXXXX)

    local log="$build_dir/cmake.log"
    if cmake -B "$build_dir" -S "$REPO_ROOT" \
        -DDFTBPLUS_ROOT="$INSTALL_PREFIX" \
        -DGRODFTB_WITH_GROMACS=OFF \
        -DCMAKE_BUILD_TYPE=Release \
        > "$log" 2>&1; then
        # Check for "Found DFTB+" status message or cached library path
        if grep -qi "Found DFTB+" "$log"; then
            pass "FindDFTBPlus.cmake found DFTB+ via DFTBPLUS_ROOT"
        elif grep -q "DFTBPlus_LIBRARY" "$build_dir/CMakeCache.txt" 2>/dev/null && \
             grep "DFTBPlus_LIBRARY" "$build_dir/CMakeCache.txt" | grep -qv "NOTFOUND"; then
            pass "DFTBPlus_LIBRARY set in CMake cache"
        else
            fail "cmake configured but DFTBPlus not found"
        fi
    else
        fail "cmake configuration failed with DFTBPLUS_ROOT"
        tail -20 "$log" 2>/dev/null
    fi

    rm -rf "$build_dir"
}

# ---------------------------------------------------------------------------
# Run all tests
# ---------------------------------------------------------------------------
echo "DFTB+ build verification tests"
echo "Install prefix: $INSTALL_PREFIX"
echo ""

test_build_script_exists
echo ""
test_dftbplus_lib_exists
echo ""
test_dftbplus_header_exists
echo ""
test_dftbplus_api_version
echo ""
test_find_dftbplus_with_root
echo ""

echo "=== Summary: $PASS passed, $FAIL failed ==="
if [[ $FAIL -gt 0 ]]; then
    exit 1
fi
exit 0
