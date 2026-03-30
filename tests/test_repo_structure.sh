#!/usr/bin/env bash
# SDD:specs.md:§4 — Verify repository directory structure matches specification.
#
# Usage: bash tests/test_repo_structure.sh [REPO_ROOT]
# Exit code: 0 if all tests pass, 1 if any fail.

set -uo pipefail

REPO_ROOT="${1:-$(cd "$(dirname "$0")/.." && pwd)}"
PASS=0
FAIL=0

pass() { PASS=$((PASS + 1)); echo "  PASS: $1"; }
fail() { FAIL=$((FAIL + 1)); echo "  FAIL: $1"; }

# ---------------------------------------------------------------------------
# test_repo_structure_dirs — AC-1: All directories from specs.md §4 exist
# ---------------------------------------------------------------------------
test_repo_structure_dirs() {
    echo "=== test_repo_structure_dirs ==="
    local dirs=(
        cmake
        external/dftbplus
        external/gromacs
        include/grodftb
        src/driver
        src/units
        src/embedding
        src/linkatom
        src/gromacs
        src/cv
        src/excited
        src/util
        tests
        tests/data
        tests/data/water_dimer
        tests/data/formamide
        tests/data/ethene
        benchmarks
        benchmarks/systems
        tools
        examples/01_gas_phase
        examples/02_solvated_qmmm
        examples/03_umbrella_et
        examples/04_surface_hopping
        docs
        containers
    )
    for d in "${dirs[@]}"; do
        if [[ -d "$REPO_ROOT/$d" ]]; then
            pass "$d exists"
        else
            fail "$d missing"
        fi
    done
}

# ---------------------------------------------------------------------------
# test_gitignore_present — AC-2: .gitignore exists and covers key patterns
# ---------------------------------------------------------------------------
test_gitignore_present() {
    echo "=== test_gitignore_present ==="
    if [[ -f "$REPO_ROOT/.gitignore" ]]; then
        pass ".gitignore exists"
    else
        fail ".gitignore missing"
        return
    fi
    if [[ -s "$REPO_ROOT/.gitignore" ]]; then
        pass ".gitignore is non-empty"
    else
        fail ".gitignore is empty"
    fi
    local patterns=("build/" "*.o" "*.so" "*.a" "*.mod" ".vscode/" ".idea/" "*.swp" ".DS_Store")
    for p in "${patterns[@]}"; do
        if grep -qF "$p" "$REPO_ROOT/.gitignore"; then
            pass ".gitignore contains '$p'"
        else
            fail ".gitignore missing pattern '$p'"
        fi
    done
}

# ---------------------------------------------------------------------------
# test_root_clean — AC-3: No unexpected files/dirs at root
# ---------------------------------------------------------------------------
test_root_clean() {
    echo "=== test_root_clean ==="
    # Allowed root-level files (from CLAUDE.md + project additions)
    local allowed_files=(
        CMakeLists.txt
        CLAUDE.md
        LICENSE
        README.md
        .gitignore
        .clang-format
        .gitmodules
        CHANGELOG.md
    )
    # Allowed root-level directories
    local allowed_dirs=(
        .git
        .claude
        cmake
        external
        include
        src
        tests
        benchmarks
        tools
        examples
        docs
        containers
    )

    # Check for unexpected files
    while IFS= read -r entry; do
        local name
        name="$(basename "$entry")"
        # Skip hidden dirs handled separately
        [[ "$name" == .git || "$name" == .claude ]] && continue
        if [[ -f "$REPO_ROOT/$name" ]]; then
            local found=false
            for a in "${allowed_files[@]}"; do
                [[ "$name" == "$a" ]] && { found=true; break; }
            done
            if $found; then
                pass "root file '$name' is allowed"
            else
                fail "unexpected root file '$name'"
            fi
        elif [[ -d "$REPO_ROOT/$name" ]]; then
            local found=false
            for a in "${allowed_dirs[@]}"; do
                [[ "$name" == "$a" ]] && { found=true; break; }
            done
            if $found; then
                pass "root dir '$name' is allowed"
            else
                fail "unexpected root dir '$name'"
            fi
        fi
    done < <(ls -A "$REPO_ROOT")
}

# ---------------------------------------------------------------------------
# test_gitkeep_present — AC-4: Empty dirs contain .gitkeep
# ---------------------------------------------------------------------------
test_gitkeep_present() {
    echo "=== test_gitkeep_present ==="
    # Directories must be tracked by git: either contain .gitkeep or real files.
    # .gitkeep is only required when the directory has no other content.
    local dirs=(
        cmake
        include/grodftb
        src/driver
        src/units
        src/embedding
        src/linkatom
        src/gromacs
        src/cv
        src/excited
        src/util
        tests/data/water_dimer
        tests/data/formamide
        tests/data/ethene
        benchmarks/systems
        tools
        examples/01_gas_phase
        examples/02_solvated_qmmm
        examples/03_umbrella_et
        examples/04_surface_hopping
        containers
    )
    for d in "${dirs[@]}"; do
        # Count non-hidden files (excluding .gitkeep itself)
        local file_count
        file_count=$(find "$REPO_ROOT/$d" -maxdepth 1 -not -name '.gitkeep' -not -name '.' -not -path "$REPO_ROOT/$d" | wc -l)
        if [[ "$file_count" -gt 0 ]]; then
            pass "$d has content (no .gitkeep needed)"
        elif [[ -f "$REPO_ROOT/$d/.gitkeep" ]]; then
            pass "$d/.gitkeep exists"
        else
            fail "$d is empty and missing .gitkeep"
        fi
    done
}

# ---------------------------------------------------------------------------
# Run all tests
# ---------------------------------------------------------------------------
echo "Repository structure tests for: $REPO_ROOT"
echo ""

test_repo_structure_dirs
echo ""
test_gitignore_present
echo ""
test_root_clean
echo ""
test_gitkeep_present
echo ""

echo "=== Summary: $PASS passed, $FAIL failed ==="
if [[ $FAIL -gt 0 ]]; then
    exit 1
fi
exit 0
