#!/usr/bin/env bash
# SDD:specs.md:§21 — Generate benchmark reference data (B1, B2, B4) by running
# standalone DFTB+ and archiving the results with full provenance.
#
# Prerequisites:
#   - DFTB+ 25.1 built at external/dftbplus/_install/bin/dftb+
#   - Slater-Koster files (mio-1-1) at tests/data/slako/mio-1-1/
#
# Usage: bash tools/generate_reference_data.sh [REPO_ROOT]

set -euo pipefail

REPO_ROOT="${1:-$(cd "$(dirname "$0")/.." && pwd)}"
DFTB_BIN="$REPO_ROOT/external/dftbplus/_install/bin/dftb+"
DFTB_LIB="$REPO_ROOT/external/dftbplus/_install/lib"
export LD_LIBRARY_PATH="${DFTB_LIB}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
SLAKO_DIR="$REPO_ROOT/tests/data/slako/mio-1-1"
DATA_DIR="$REPO_ROOT/tests/data"
TESTPARAMS_COMMIT="0e1d95abc70339d017fd53d40c41c041057fd03f"
TESTPARAMS_URL="https://github.com/dftbplus/testparams/archive/${TESTPARAMS_COMMIT}.tar.gz"

# ---------------------------------------------------------------------------
# Preflight checks
# ---------------------------------------------------------------------------
if [[ ! -x "$DFTB_BIN" ]]; then
    echo "ERROR: DFTB+ binary not found at $DFTB_BIN"
    echo "       Build DFTB+ first: bash tools/build_dftbplus.sh"
    exit 1
fi

# ---------------------------------------------------------------------------
# Step 0: Download Slater-Koster files if not present
# ---------------------------------------------------------------------------
download_slako() {
    if [[ -f "$SLAKO_DIR/O-O.skf" ]]; then
        echo "SK files already present at $SLAKO_DIR"
        return 0
    fi
    echo "Downloading mio-1-1 Slater-Koster parameter files..."
    local tmpdir
    tmpdir=$(mktemp -d)
    local tarball="$tmpdir/testparams.tar.gz"

    if ! curl -fsSL -o "$tarball" "$TESTPARAMS_URL"; then
        echo "ERROR: Failed to download test parameters from $TESTPARAMS_URL"
        echo "       You may need to download them manually."
        echo "       Expected location: $SLAKO_DIR/"
        echo "       Source: https://github.com/dftbplus/testparams"
        rm -rf "$tmpdir"
        exit 1
    fi

    tar xzf "$tarball" -C "$tmpdir"
    local extracted="$tmpdir/testparams-${TESTPARAMS_COMMIT}"

    mkdir -p "$SLAKO_DIR"
    cp "$extracted/slakos/origin/mio-1-1/"*.skf "$SLAKO_DIR/"

    rm -rf "$tmpdir"
    echo "SK files installed at $SLAKO_DIR"
}

# ---------------------------------------------------------------------------
# Helper: get DFTB+ version string
# ---------------------------------------------------------------------------
get_dftb_version() {
    "$DFTB_BIN" --version 2>/dev/null | grep -oP 'DFTB\+.*' | head -1 || echo "unknown"
}

# ---------------------------------------------------------------------------
# Helper: extract DFTB+ git commit from version output
# ---------------------------------------------------------------------------
get_dftb_commit() {
    local version_out
    version_out=$("$DFTB_BIN" --version 2>/dev/null || true)
    # Try to extract commit hash from version output
    echo "$version_out" | grep -oE '[0-9a-f]{7,40}' | head -1 || echo "unknown"
}

# ---------------------------------------------------------------------------
# Helper: parse autotest.tag for energy, forces, charges
# ---------------------------------------------------------------------------
parse_autotest_tag() {
    local tag_file="$1"
    local json_file="$2"
    local natom="$3"

    if [[ ! -f "$tag_file" ]]; then
        echo "ERROR: autotest.tag not found at $tag_file"
        return 1
    fi

    python3 - "$tag_file" "$json_file" "$natom" <<'PYEOF'
import sys
import json
import re

tag_file = sys.argv[1]
json_file = sys.argv[2]
natom = int(sys.argv[3])

with open(tag_file) as f:
    content = f.read()

result = {}

# Parse total energy (Mermin free energy)
# Format: mermin_energy        :real:0:
#   -4.0709442785266400E+00
m = re.search(r'mermin_energy\s*:real:0:\s*\n\s*([^\n]+)', content)
if m:
    result["energy_hartree"] = float(m.group(1).strip())

# Also try total_energy if mermin not found
if "energy_hartree" not in result:
    m = re.search(r'total_energy\s*:real:0:\s*\n\s*([^\n]+)', content)
    if m:
        result["energy_hartree"] = float(m.group(1).strip())

# Parse forces
# Format: forces               :real:2:3,N
#   val1 val2 val3 ...
m = re.search(r'forces\s*:real:2:3,' + str(natom) + r'\s*\n((?:\s*[^\n]+\n?)+?)(?=\w|\Z)', content)
if m:
    force_text = m.group(1).strip()
    force_vals = [float(x) for x in force_text.split()]
    # Reshape into [natom][3]
    forces = []
    for i in range(natom):
        forces.append(force_vals[i*3:(i+1)*3])
    result["forces_hartree_bohr"] = forces

# Parse gross charges (try gross_charges first, then sum orbital_charges)
# Format: gross_charges         :real:1:N
m = re.search(r'gross_charges\s*:real:1:' + str(natom) + r'\s*\n((?:\s*[^\n]+\n?)+?)(?=\w|\Z)', content)
if m:
    charge_text = m.group(1).strip()
    result["gross_charges_e"] = [float(x) for x in charge_text.split()]
else:
    # Try parsing from detailed.out (net Mulliken charges)
    import os
    det_file = os.path.join(os.path.dirname(tag_file), "detailed.out")
    if os.path.exists(det_file):
        with open(det_file) as df:
            det = df.read()
        m2 = re.search(r'Atomic gross charges \(e\)\s*\n\s*Atom\s+Charge\s*\n((?:\s+\d+\s+[^\n]+\n)+)', det)
        if m2:
            lines = m2.group(1).strip().split('\n')
            charges = [float(l.split()[1]) for l in lines]
            result["gross_charges_e"] = charges

with open(json_file, 'w') as f:
    json.dump(result, f, indent=2)

print(f"Wrote reference data to {json_file}")
if "energy_hartree" in result:
    print(f"  Energy: {result['energy_hartree']:.16e} Ha")
if "forces_hartree_bohr" in result:
    print(f"  Forces: {len(result['forces_hartree_bohr'])} atoms")
if "gross_charges_e" in result:
    print(f"  Charges: {len(result['gross_charges_e'])} atoms")
PYEOF
}

# ---------------------------------------------------------------------------
# Helper: write provenance file
# ---------------------------------------------------------------------------
write_provenance() {
    local dest="$1"
    local benchmark="$2"
    local dftb_version
    dftb_version=$(get_dftb_version)
    local dftb_commit
    dftb_commit=$(get_dftb_commit)

    cat > "$dest" <<EOF
Benchmark: $benchmark
Date: $(date -u +"%Y-%m-%dT%H:%M:%SZ")
DFTB+ version: $dftb_version
DFTB+ commit: $dftb_commit
DFTB+ binary: $DFTB_BIN
SK parameter set: mio-1-1
SK source: https://github.com/dftbplus/testparams (commit $TESTPARAMS_COMMIT)
Generator script: tools/generate_reference_data.sh
Host: $(hostname)
User: $(whoami)
EOF
    echo "Wrote provenance to $dest"
}

# ---------------------------------------------------------------------------
# Run a benchmark
# ---------------------------------------------------------------------------
run_benchmark() {
    local name="$1"
    local natom="$2"
    local bench_dir="$DATA_DIR/$name"

    echo ""
    echo "=== Running benchmark $name ==="

    if [[ ! -f "$bench_dir/dftb_in.hsd" ]]; then
        echo "ERROR: $bench_dir/dftb_in.hsd not found"
        return 1
    fi
    if [[ ! -f "$bench_dir/geo.gen" ]]; then
        echo "ERROR: $bench_dir/geo.gen not found"
        return 1
    fi

    # Run DFTB+ in the benchmark directory
    local workdir
    workdir=$(mktemp -d)
    cp "$bench_dir/dftb_in.hsd" "$workdir/"
    cp "$bench_dir/geo.gen" "$workdir/"

    # Create symlink to SK files (the HSD references ../slako/mio-1-1/)
    mkdir -p "$workdir/../slako"
    ln -sf "$SLAKO_DIR" "$workdir/../slako/mio-1-1"

    echo "Running DFTB+ in $workdir ..."
    local dftb_output
    if ! dftb_output=$(cd "$workdir" && "$DFTB_BIN" 2>&1); then
        echo "ERROR: DFTB+ failed for $name"
        echo "$dftb_output"
        rm -rf "$workdir"
        return 1
    fi
    echo "DFTB+ completed successfully"

    # Save the full DFTB+ output
    echo "$dftb_output" > "$bench_dir/dftb_output.log"

    # Copy autotest.tag
    if [[ -f "$workdir/autotest.tag" ]]; then
        cp "$workdir/autotest.tag" "$bench_dir/autotest.tag"
    fi

    # Copy detailed.out if present
    if [[ -f "$workdir/detailed.out" ]]; then
        cp "$workdir/detailed.out" "$bench_dir/detailed.out"
    fi

    # Parse results into reference.json
    parse_autotest_tag "$workdir/autotest.tag" "$bench_dir/reference.json" "$natom"

    # Write provenance
    write_provenance "$bench_dir/provenance.txt" "$name"

    rm -rf "$workdir"
    echo "=== $name complete ==="
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
echo "GRO-DFTB Reference Data Generator"
echo "================================="
echo "DFTB+ binary: $DFTB_BIN"
echo ""

download_slako

run_benchmark "b1" 6    # Water dimer
run_benchmark "b2" 6    # Formamide
run_benchmark "b4" 6    # Water dimer + point charge

echo ""
echo "All benchmarks complete."
echo "Reference data written to $DATA_DIR/b{1,2,4}/reference.json"
echo "Provenance files written to $DATA_DIR/b{1,2,4}/provenance.txt"
