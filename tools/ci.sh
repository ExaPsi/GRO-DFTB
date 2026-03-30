#!/usr/bin/env bash
# Local CI pipeline for GRO-DFTB (solo developer workflow)
# Usage: tools/ci.sh [--clean] [--sanitize] [--verbose]
#
# Runs: configure → build → test → (optional) sanitizer pass
# Returns: 0 on success, 1 on failure

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD_DIR="${REPO_ROOT}/_build"

# Parse flags
CLEAN=false
SANITIZE=false
VERBOSE=false
for arg in "$@"; do
    case "$arg" in
        --clean)    CLEAN=true ;;
        --sanitize) SANITIZE=true ;;
        --verbose)  VERBOSE=true ;;
        --help|-h)
            echo "Usage: tools/ci.sh [--clean] [--sanitize] [--verbose]"
            echo "  --clean     Remove build directory before building"
            echo "  --sanitize  Additional build+test pass with ASan+UBSan"
            echo "  --verbose   Show full build and test output"
            exit 0
            ;;
        *)
            echo "Unknown flag: $arg" >&2
            exit 1
            ;;
    esac
done

# Counters
STEPS_RUN=0
STEPS_PASSED=0
STEPS_FAILED=0

step_pass() {
    STEPS_RUN=$((STEPS_RUN + 1))
    STEPS_PASSED=$((STEPS_PASSED + 1))
    echo "  PASS: $1"
}

step_fail() {
    STEPS_RUN=$((STEPS_RUN + 1))
    STEPS_FAILED=$((STEPS_FAILED + 1))
    echo "  FAIL: $1" >&2
}

echo "=== GRO-DFTB Local CI ==="
echo "Repository: ${REPO_ROOT}"
echo "Build dir:  ${BUILD_DIR}"
echo ""

# --- Clean ---
if $CLEAN && [ -d "$BUILD_DIR" ]; then
    echo "[1/4] Cleaning build directory..."
    rm -rf "$BUILD_DIR"
fi

# --- Configure ---
echo "[1/4] Configuring..."
# Auto-detect DFTB+ install if available
DFTBPLUS_INSTALL="${REPO_ROOT}/external/dftbplus/_install"
DFTBPLUS_FLAG=""
if [ -d "$DFTBPLUS_INSTALL" ]; then
    DFTBPLUS_FLAG="-DDFTBPLUS_ROOT=$DFTBPLUS_INSTALL"
fi

CMAKE_FLAGS=(
    -B "$BUILD_DIR"
    -S "$REPO_ROOT"
    -DGRODFTB_WITH_GROMACS=OFF
    -DCMAKE_C_FLAGS="-Wall -Wextra -pedantic -Werror"
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
    ${DFTBPLUS_FLAG:+$DFTBPLUS_FLAG}
)

CONFIG_RC=0
if $VERBOSE; then
    cmake "${CMAKE_FLAGS[@]}" 2>&1 || CONFIG_RC=$?
else
    cmake "${CMAKE_FLAGS[@]}" > /dev/null 2>&1 || CONFIG_RC=$?
fi

if [ $CONFIG_RC -eq 0 ]; then
    step_pass "CMake configure"
else
    step_fail "CMake configure"
    echo ""
    echo "=== CI SUMMARY: ${STEPS_PASSED}/${STEPS_RUN} passed, ${STEPS_FAILED} failed — FAIL ==="
    exit 1
fi

# --- Build ---
echo "[2/4] Building..."
BUILD_RC=0
if $VERBOSE; then
    cmake --build "$BUILD_DIR" 2>&1 || BUILD_RC=$?
else
    BUILD_OUTPUT=$(cmake --build "$BUILD_DIR" 2>&1) || BUILD_RC=$?
fi

if [ $BUILD_RC -eq 0 ]; then
    step_pass "Build (zero warnings with -Werror)"
else
    step_fail "Build"
    if ! $VERBOSE; then
        echo "$BUILD_OUTPUT" >&2
    fi
    echo ""
    echo "=== CI SUMMARY: ${STEPS_PASSED}/${STEPS_RUN} passed, ${STEPS_FAILED} failed — FAIL ==="
    exit 1
fi

# --- Test ---
echo "[3/4] Testing..."
TEST_RC=0
if $VERBOSE; then
    ctest --test-dir "$BUILD_DIR" --output-on-failure 2>&1 || TEST_RC=$?
else
    TEST_OUTPUT=$(ctest --test-dir "$BUILD_DIR" --output-on-failure 2>&1) || TEST_RC=$?
fi

if [ $TEST_RC -eq 0 ]; then
    # Extract test counts from ctest output
    if $VERBOSE; then
        step_pass "Tests"
    else
        TCOUNT=$(echo "$TEST_OUTPUT" | grep -oP '\d+ tests? passed' || echo "tests passed")
        step_pass "Tests ($TCOUNT)"
    fi
else
    step_fail "Tests"
    if ! $VERBOSE; then
        echo "$TEST_OUTPUT" >&2
    fi
    echo ""
    echo "=== CI SUMMARY: ${STEPS_PASSED}/${STEPS_RUN} passed, ${STEPS_FAILED} failed — FAIL ==="
    exit 1
fi

# --- Sanitizer pass (optional) ---
if $SANITIZE; then
    echo "[4/4] Sanitizer pass (ASan + UBSan)..."
    SANBUILD="${BUILD_DIR}_sanitize"
    rm -rf "$SANBUILD"

    SAN_FLAGS=(
        -B "$SANBUILD"
        -S "$REPO_ROOT"
        -DGRODFTB_WITH_GROMACS=OFF
        -DCMAKE_C_FLAGS="-Wall -Wextra -pedantic -fsanitize=address,undefined -fno-omit-frame-pointer"
        -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined"
        -DCMAKE_SHARED_LINKER_FLAGS="-fsanitize=address,undefined"
        ${DFTBPLUS_FLAG:+$DFTBPLUS_FLAG}
    )

    SAN_RC=0
    if $VERBOSE; then
        { cmake "${SAN_FLAGS[@]}" 2>&1 && \
          cmake --build "$SANBUILD" 2>&1 && \
          ctest --test-dir "$SANBUILD" --output-on-failure 2>&1; } || SAN_RC=$?
    else
        SAN_OUTPUT=$({ cmake "${SAN_FLAGS[@]}" 2>&1 && \
                       cmake --build "$SANBUILD" 2>&1 && \
                       ctest --test-dir "$SANBUILD" --output-on-failure 2>&1; }) || SAN_RC=$?
    fi

    if [ $SAN_RC -eq 0 ]; then
        step_pass "Sanitizers (ASan + UBSan)"
    else
        step_fail "Sanitizers (ASan + UBSan)"
        if ! $VERBOSE && [ -n "${SAN_OUTPUT:-}" ]; then
            echo "$SAN_OUTPUT" >&2
        fi
    fi

    rm -rf "$SANBUILD"
else
    echo "[4/4] Sanitizer pass — skipped (use --sanitize to enable)"
fi

# --- Summary ---
echo ""
if [ $STEPS_FAILED -eq 0 ]; then
    echo "=== CI SUMMARY: ${STEPS_PASSED}/${STEPS_RUN} passed, 0 failed — PASS ==="
    exit 0
else
    echo "=== CI SUMMARY: ${STEPS_PASSED}/${STEPS_RUN} passed, ${STEPS_FAILED} failed — FAIL ==="
    exit 1
fi
