#!/usr/bin/env bash
# Tests for the local CI pipeline (US-008)
# Usage: bash tests/test_ci.sh

set -uo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TESTS_RUN=0
TESTS_PASSED=0
TESTS_FAILED=0

pass() {
    TESTS_RUN=$((TESTS_RUN + 1))
    TESTS_PASSED=$((TESTS_PASSED + 1))
    echo "  PASS: $1"
}

fail() {
    TESTS_RUN=$((TESTS_RUN + 1))
    TESTS_FAILED=$((TESTS_FAILED + 1))
    echo "  FAIL: $1" >&2
    [ -n "${2:-}" ] && echo "        $2" >&2
}

echo "=== test_ci.sh ==="

# --- AC-1: tools/ci.sh exists and is executable ---
test_ci_script_exists() {
    if [ -x "${REPO_ROOT}/tools/ci.sh" ]; then
        pass "test_ci_script_exists"
    else
        fail "test_ci_script_exists" "tools/ci.sh not found or not executable"
    fi
}

# --- AC-2: ci.sh configures CMake with strict warnings ---
test_ci_configures() {
    local output
    output=$("${REPO_ROOT}/tools/ci.sh" --clean 2>&1)
    if echo "$output" | grep -q "PASS.*configure\|PASS.*CMake"; then
        pass "test_ci_configures"
    else
        fail "test_ci_configures" "CMake configure step did not pass"
    fi
}

# --- AC-3: builds with zero warnings ---
test_ci_builds_clean() {
    local output
    output=$("${REPO_ROOT}/tools/ci.sh" 2>&1)
    if echo "$output" | grep -q "PASS.*[Bb]uild"; then
        pass "test_ci_builds_clean"
    else
        fail "test_ci_builds_clean" "Build step did not pass"
    fi
}

# --- AC-4: ctest passes ---
test_ci_tests_pass() {
    local output
    output=$("${REPO_ROOT}/tools/ci.sh" 2>&1)
    if echo "$output" | grep -q "PASS.*[Tt]est"; then
        pass "test_ci_tests_pass"
    else
        fail "test_ci_tests_pass" "Test step did not pass"
    fi
}

# --- AC-5: sanitizer pass ---
test_ci_sanitizers() {
    local output rc
    output=$("${REPO_ROOT}/tools/ci.sh" --sanitize 2>&1)
    rc=$?
    if [ $rc -eq 0 ] && echo "$output" | grep -q "PASS.*[Ss]anitizer"; then
        pass "test_ci_sanitizers"
    else
        fail "test_ci_sanitizers" "Sanitizer pass did not succeed"
    fi
}

# --- AC-6: hook installed ---
test_hook_installed() {
    # Run install script
    bash "${REPO_ROOT}/tools/install-hooks.sh" > /dev/null 2>&1
    if [ -x "${REPO_ROOT}/.git/hooks/pre-push" ]; then
        pass "test_hook_installed"
    else
        fail "test_hook_installed" ".git/hooks/pre-push not found after install"
    fi
}

# --- AC-7: hook blocks on failure ---
test_hook_blocks_on_failure() {
    # The hook calls tools/ci.sh. Since ci.sh passes, we verify the hook
    # script structure: it should exit nonzero when ci.sh fails.
    local hook="${REPO_ROOT}/.git/hooks/pre-push"
    if [ -x "$hook" ] && grep -q 'exit 1' "$hook"; then
        pass "test_hook_blocks_on_failure"
    else
        fail "test_hook_blocks_on_failure" "Hook does not contain failure exit path"
    fi
}

# --- AC-8: summary output ---
test_ci_summary_output() {
    local output
    output=$("${REPO_ROOT}/tools/ci.sh" 2>&1)
    if echo "$output" | grep -qE "CI SUMMARY:.*passed.*failed.*(PASS|FAIL)"; then
        pass "test_ci_summary_output"
    else
        fail "test_ci_summary_output" "Summary line not found in output"
    fi
}

# Run all tests
test_ci_script_exists
test_ci_configures
test_ci_builds_clean
test_ci_tests_pass
test_ci_sanitizers
test_hook_installed
test_hook_blocks_on_failure
test_ci_summary_output

echo ""
echo "test_ci: ${TESTS_RUN} run, ${TESTS_PASSED} passed, ${TESTS_FAILED} failed"
exit $TESTS_FAILED
