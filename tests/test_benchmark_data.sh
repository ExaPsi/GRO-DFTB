#!/usr/bin/env bash
# SDD:specs.md:§21 — Verify benchmark reference data (B1, B2, B4) exists,
# is complete, has provenance, and is reproducible.
#
# Usage: bash tests/test_benchmark_data.sh [REPO_ROOT]
# Exit code: 0 if all tests pass, 1 if any fail.

set -uo pipefail

REPO_ROOT="${1:-$(cd "$(dirname "$0")/.." && pwd)}"
DATA_DIR="$REPO_ROOT/tests/data"
DFTB_BIN="$REPO_ROOT/external/dftbplus/_install/bin/dftb+"
SLAKO_DIR="$REPO_ROOT/tests/data/slako/mio-1-1"
PASS=0
FAIL=0
SKIP=0

pass() { PASS=$((PASS + 1)); echo "  PASS: $1"; }
fail() { FAIL=$((FAIL + 1)); echo "  FAIL: $1"; }
skip() { SKIP=$((SKIP + 1)); echo "  SKIP: $1"; }

# ---------------------------------------------------------------------------
# test_b1_input_files_exist — B1 water dimer input files present
# ---------------------------------------------------------------------------
test_b1_input_files_exist() {
    echo "=== test_b1_input_files_exist ==="
    local ok=1
    if [[ ! -f "$DATA_DIR/b1/dftb_in.hsd" ]]; then
        fail "tests/data/b1/dftb_in.hsd missing"
        ok=0
    fi
    if [[ ! -f "$DATA_DIR/b1/geo.gen" ]]; then
        fail "tests/data/b1/geo.gen missing"
        ok=0
    fi
    if [[ $ok -eq 1 ]]; then
        pass "B1 input files (dftb_in.hsd, geo.gen) present"
    fi
}

# ---------------------------------------------------------------------------
# test_b2_input_files_exist — B2 formamide input files present
# ---------------------------------------------------------------------------
test_b2_input_files_exist() {
    echo "=== test_b2_input_files_exist ==="
    local ok=1
    if [[ ! -f "$DATA_DIR/b2/dftb_in.hsd" ]]; then
        fail "tests/data/b2/dftb_in.hsd missing"
        ok=0
    fi
    if [[ ! -f "$DATA_DIR/b2/geo.gen" ]]; then
        fail "tests/data/b2/geo.gen missing"
        ok=0
    fi
    if [[ $ok -eq 1 ]]; then
        pass "B2 input files (dftb_in.hsd, geo.gen) present"
    fi
}

# ---------------------------------------------------------------------------
# test_b4_input_files_exist — B4 water dimer + point charge input files present
# ---------------------------------------------------------------------------
test_b4_input_files_exist() {
    echo "=== test_b4_input_files_exist ==="
    local ok=1
    if [[ ! -f "$DATA_DIR/b4/dftb_in.hsd" ]]; then
        fail "tests/data/b4/dftb_in.hsd missing"
        ok=0
    fi
    if [[ ! -f "$DATA_DIR/b4/geo.gen" ]]; then
        fail "tests/data/b4/geo.gen missing"
        ok=0
    fi
    if [[ $ok -eq 1 ]]; then
        pass "B4 input files (dftb_in.hsd, geo.gen) present"
    fi
}

# ---------------------------------------------------------------------------
# Helper: check that reference.json has required fields
# ---------------------------------------------------------------------------
check_reference_json() {
    local json_file="$1"
    local label="$2"

    if [[ ! -f "$json_file" ]]; then
        fail "$label: reference.json missing"
        return 1
    fi

    local ok=1
    for field in energy_hartree forces_hartree_bohr gross_charges_e; do
        if ! python3 -c "
import json, sys
with open('$json_file') as f:
    d = json.load(f)
if '$field' not in d:
    sys.exit(1)
if d['$field'] is None:
    sys.exit(1)
" 2>/dev/null; then
            fail "$label: reference.json missing field '$field'"
            ok=0
        fi
    done

    if [[ $ok -eq 1 ]]; then
        pass "$label: reference.json has energy, forces, charges"
    fi
}

# ---------------------------------------------------------------------------
# test_b1_reference_data_exists — B1 reference.json has required fields
# ---------------------------------------------------------------------------
test_b1_reference_data_exists() {
    echo "=== test_b1_reference_data_exists ==="
    check_reference_json "$DATA_DIR/b1/reference.json" "B1"
}

# ---------------------------------------------------------------------------
# test_b2_reference_data_exists — B2 reference.json has required fields
# ---------------------------------------------------------------------------
test_b2_reference_data_exists() {
    echo "=== test_b2_reference_data_exists ==="
    check_reference_json "$DATA_DIR/b2/reference.json" "B2"
}

# ---------------------------------------------------------------------------
# test_b4_reference_data_exists — B4 reference.json has required fields
# ---------------------------------------------------------------------------
test_b4_reference_data_exists() {
    echo "=== test_b4_reference_data_exists ==="
    check_reference_json "$DATA_DIR/b4/reference.json" "B4"
}

# ---------------------------------------------------------------------------
# test_provenance_complete — each benchmark has provenance with required fields
# ---------------------------------------------------------------------------
test_provenance_complete() {
    echo "=== test_provenance_complete ==="
    local ok=1
    for bm in b1 b2 b4; do
        local prov="$DATA_DIR/$bm/provenance.txt"
        if [[ ! -f "$prov" ]]; then
            fail "$bm: provenance.txt missing"
            ok=0
            continue
        fi
        for field in "version" "commit" "date"; do
            if ! grep -qi "$field" "$prov"; then
                fail "$bm: provenance.txt missing '$field' field"
                ok=0
            fi
        done
    done
    if [[ $ok -eq 1 ]]; then
        pass "All provenance files complete (version, commit, date)"
    fi
}

# ---------------------------------------------------------------------------
# test_dftbplus_reproducible — re-run DFTB+ on B1 and compare energy
# ---------------------------------------------------------------------------
test_dftbplus_reproducible() {
    echo "=== test_dftbplus_reproducible ==="

    # Skip if DFTB+ binary or SK files not available
    if [[ ! -x "$DFTB_BIN" ]]; then
        skip "DFTB+ binary not found at $DFTB_BIN"
        return 0
    fi
    if [[ ! -f "$SLAKO_DIR/O-O.skf" ]]; then
        skip "SK files not found at $SLAKO_DIR"
        return 0
    fi
    if [[ ! -f "$DATA_DIR/b1/reference.json" ]]; then
        skip "B1 reference.json not found — run tools/generate_reference_data.sh first"
        return 0
    fi

    # Get archived energy
    local archived_energy
    archived_energy=$(python3 -c "
import json
with open('$DATA_DIR/b1/reference.json') as f:
    d = json.load(f)
print(f\"{d['energy_hartree']:.16e}\")
" 2>/dev/null)

    if [[ -z "$archived_energy" ]]; then
        fail "Could not read archived energy from B1 reference.json"
        return 1
    fi

    # Run DFTB+ in a temporary directory
    local workdir
    workdir=$(mktemp -d)
    cp "$DATA_DIR/b1/dftb_in.hsd" "$workdir/"
    cp "$DATA_DIR/b1/geo.gen" "$workdir/"
    mkdir -p "$workdir/../slako"
    ln -sf "$SLAKO_DIR" "$workdir/../slako/mio-1-1"

    local dftb_out
    if ! dftb_out=$(cd "$workdir" && "$DFTB_BIN" 2>&1); then
        fail "DFTB+ failed during reproducibility check"
        echo "  Output: $dftb_out"
        rm -rf "$workdir"
        return 1
    fi

    if [[ ! -f "$workdir/autotest.tag" ]]; then
        fail "autotest.tag not generated during reproducibility check"
        rm -rf "$workdir"
        return 1
    fi

    # Extract energy from fresh run
    local fresh_energy
    fresh_energy=$(python3 -c "
import re
with open('$workdir/autotest.tag') as f:
    content = f.read()
m = re.search(r'mermin_energy\s*:real:0:\s*\n\s*([^\n]+)', content)
if m:
    print(f'{float(m.group(1).strip()):.16e}')
else:
    m = re.search(r'total_energy\s*:real:0:\s*\n\s*([^\n]+)', content)
    if m:
        print(f'{float(m.group(1).strip()):.16e}')
" 2>/dev/null)

    rm -rf "$workdir"

    if [[ -z "$fresh_energy" ]]; then
        fail "Could not extract energy from fresh DFTB+ run"
        return 1
    fi

    # Compare energies — must match to machine precision (< 1e-14 Ha)
    local match
    match=$(python3 -c "
import sys
a = float('$archived_energy')
b = float('$fresh_energy')
diff = abs(a - b)
if diff < 1.0e-14:
    print('MATCH')
    print(f'  Archived: {a:.16e} Ha')
    print(f'  Fresh:    {b:.16e} Ha')
    print(f'  Diff:     {diff:.2e} Ha')
else:
    print('MISMATCH')
    print(f'  Archived: {a:.16e} Ha')
    print(f'  Fresh:    {b:.16e} Ha')
    print(f'  Diff:     {diff:.2e} Ha')
" 2>/dev/null)

    if echo "$match" | grep -q "MATCH"; then
        pass "B1 energy reproducible to machine precision"
        echo "$match" | tail -3 | while read -r line; do echo "    $line"; done
    else
        fail "B1 energy NOT reproducible"
        echo "$match" | tail -3 | while read -r line; do echo "    $line"; done
    fi
}

# ---------------------------------------------------------------------------
# Run all tests
# ---------------------------------------------------------------------------
echo "Benchmark reference data verification tests"
echo "Data directory: $DATA_DIR"
echo ""

test_b1_input_files_exist
echo ""
test_b2_input_files_exist
echo ""
test_b4_input_files_exist
echo ""
test_b1_reference_data_exists
echo ""
test_b2_reference_data_exists
echo ""
test_b4_reference_data_exists
echo ""
test_provenance_complete
echo ""
test_dftbplus_reproducible
echo ""

echo "=== Summary: $PASS passed, $FAIL failed, $SKIP skipped ==="
if [[ $FAIL -gt 0 ]]; then
    exit 1
fi
exit 0
