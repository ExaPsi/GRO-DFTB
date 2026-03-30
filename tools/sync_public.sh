#!/bin/bash
# sync_public.sh — Sync minimal necessary files to public release directory
#
# Usage: bash tools/sync_public.sh
#
# Creates/updates -public/ with only the files
# needed for the public GitHub release at https://github.com/ExaPsi/GRO-DFTB
#
# Re-runnable: safe to run multiple times; --delete removes stale files.

set -euo pipefail

SRC=""
DST="-public"

echo "=== GRO-DFTB Public Release Sync ==="
echo "Source: $SRC"
echo "Destination: $DST"
echo ""

# ── Step 1: rsync with exclusions ────────────────────────────────────────

mkdir -p "$DST"

rsync -av --delete \
  --exclude='.git/' \
  --exclude='.claude/' \
  --exclude='CLAUDE.md' \
  --exclude='docs/' \
  --exclude='build/' \
  --exclude='_build/' \
  --exclude='external/' \
  --exclude='*.pyc' \
  --exclude='__pycache__/' \
  --exclude='examples/01_gas_phase/' \
  --exclude='examples/02_solvated_qmmm/' \
  --exclude='examples/03_umbrella_et/' \
  --exclude='examples/04_surface_hopping/' \
  --exclude='benchmarks/systems/' \
  --exclude='containers/' \
  --exclude='tests/data/b5/*.xvg' \
  --exclude='tests/data/b5/old_run_feb4/' \
  --exclude='tests/data/b5/drift_analysis*.json' \
  --exclude='tests/data/b5/provenance_template.json' \
  --exclude='tests/data/b5/dftb_template.hsd' \
  --exclude='tests/data/b5/dftb_mio_template.hsd' \
  --exclude='tests/data/b5/nve_damped/attempt1_failed/' \
  --exclude='tests/data/b5/nve_damped/*.edr' \
  --exclude='tests/data/b5/nve_damped/*.xtc' \
  --exclude='tests/data/b5/nve_damped/*.trr' \
  --exclude='tests/data/b5/nve_damped/*.cpt' \
  --exclude='tests/data/b5/nve_damped/*.tpr' \
  --exclude='tests/data/b5/nve_damped/*.log' \
  --exclude='tests/data/b5/nve_damped/#*' \
  --exclude='tests/data/b5/*.edr' \
  --exclude='tests/data/b5/*.xtc' \
  --exclude='tests/data/b5/*.trr' \
  --exclude='tests/data/b5/*.cpt' \
  --exclude='tests/data/b5/*.tpr' \
  --exclude='tests/data/b5/#*' \
  --exclude='tests/data/b5/b5_nve_500ps*' \
  --exclude='tests/data/b5/mdrun_500ps*' \
  --exclude='tests/data/b6_solvated/*.xvg' \
  --exclude='tests/data/b5_energy_decomposition/*.xvg' \
  --exclude='tests/data/b4_damped/*.bak' \
  --exclude='tests/data/b4/fresh_output.log' \
  --exclude='tests/data/water_dimer/' \
  --exclude='tests/data/formamide/' \
  --exclude='tests/data/ethene/' \
  --exclude='#*#' \
  --exclude='*.swp' \
  --exclude='*.swo' \
  --exclude='*~' \
  "$SRC/" "$DST/"

echo ""
echo "=== Step 1 complete: files synced ==="

# ── Step 2: Remove circular symlinks ─────────────────────────────────────

for link in "$DST/tests/data/slako/slako" \
            "$DST/tests/data/slako/mio-1-1/mio-1-1"; do
  if [ -L "$link" ]; then
    rm -f "$link"
    echo "Removed circular symlink: $link"
  fi
done

echo ""
echo "=== Step 2 complete: symlinks cleaned ==="

# ── Step 3: Fix hardcoded paths ──────────────────────────────────────────

echo "Fixing hardcoded paths..."

# Replace  with relative or generic paths
find "$DST" -type f \( \
  -name "*.hsd" -o -name "*.mdp" -o -name "*.json" -o \
  -name "*.txt" -o -name "*.py" -o -name "*.c" -o \
  -name "*.md" -o -name "*.sh" \
\) -exec sed -i \
  -e 's|gmx|gmx|g' \
  -e 's|dftb+|dftb+|g' \
  -e 's|../slako/|../slako/|g' \
  -e 's|||g' \
  -e 's|||g' \
  {} +

# Also fix in .xvg file headers (GROMACS comment lines)
find "$DST" -type f -name "*.xvg" -exec sed -i \
  -e 's|[^ ]*|<build-dir>|g' \
  -e 's|/home/vyv/[^ ]*|<local>|g' \
  {} +

echo ""
echo "=== Step 3 complete: paths fixed ==="

# ── Step 4: Update GitHub URL ────────────────────────────────────────────

echo "Updating GitHub URLs..."

find "$DST" -type f \( -name "*.md" -o -name "*.tex" -o -name "*.cff" \) \
  -exec sed -i 's|https://github.com/yyods/GRO-DFTB|https://github.com/ExaPsi/GRO-DFTB|g' {} +

echo ""
echo "=== Step 4 complete: URLs updated ==="

# ── Step 5: Create CITATION.cff ──────────────────────────────────────────

cat > "$DST/CITATION.cff" << 'CITEOF'
cff-version: 1.2.0
message: "If you use this software, please cite the article below."
title: "GRO-DFTB: Integrating SCC-DFTB with GROMACS for Energy-Conserving QM/MM Molecular Dynamics"
type: software
authors:
  - family-names: Vchirawongkwin
    given-names: Viwat
    orcid: "https://orcid.org/0000-0002-2747-1552"
    affiliation: "Chulalongkorn University"
version: "1.0"
license: LGPL-3.0-or-later
repository-code: "https://github.com/ExaPsi/GRO-DFTB"
preferred-citation:
  type: article
  title: "GRO-DFTB: Integrating SCC-DFTB with GROMACS for Energy-Conserving QM/MM Molecular Dynamics"
  authors:
    - family-names: Vchirawongkwin
      given-names: Viwat
  journal: "Journal of Chemical Theory and Computation"
  year: 2027
CITEOF

echo "Created CITATION.cff"

# ── Step 6: Create CONTRIBUTING.md ───────────────────────────────────────

cat > "$DST/CONTRIBUTING.md" << 'CONTEOF'
# Contributing to GRO-DFTB

Thank you for your interest in contributing to GRO-DFTB.

## Reporting Issues

Please report bugs and feature requests via [GitHub Issues](https://github.com/ExaPsi/GRO-DFTB/issues).

When reporting a bug, include:
- GRO-DFTB version (`git describe --tags`)
- DFTB+ and GROMACS versions
- Operating system and compiler
- Minimal reproducing input files
- Full error output

## Contributing Code

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-feature`)
3. Write tests for your changes
4. Ensure all existing tests pass (`ctest --test-dir build`)
5. Submit a pull request

### Code Style

- C11 standard for `libgrodftb` (src/, include/)
- C++17 for GROMACS integration code
- Function names: `grodftb_` prefix for public API
- All public functions documented in headers

### Testing Requirements

- Every new function must have a corresponding test
- Finite-difference validation for any force-producing code path
- NVE energy conservation test for new embedding modes

## License

All contributions are licensed under LGPL-3.0-or-later, consistent with the project license.
CONTEOF

echo "Created CONTRIBUTING.md"

# ── Step 7: Report ───────────────────────────────────────────────────────

echo ""
echo "=== Sync complete ==="
echo ""
echo "File count:"
find "$DST" -type f | wc -l
echo ""
echo "Directory size:"
du -sh "$DST"
echo ""
echo "Size breakdown:"
du -sh "$DST/src" "$DST/include" "$DST/tests" "$DST/tools" "$DST/examples" "$DST/cmake" 2>/dev/null
echo ""

# Verify no private paths remain
echo "Checking for remaining /home/vyv/ paths..."
FOUND=$(grep -rl "/home/vyv/" "$DST" --include="*.c" --include="*.h" --include="*.hsd" --include="*.mdp" --include="*.json" --include="*.py" --include="*.txt" --include="*.md" 2>/dev/null || true)
if [ -z "$FOUND" ]; then
  echo "  PASS: No /home/vyv/ paths found"
else
  echo "  FAIL: Found paths in:"
  echo "$FOUND"
fi

echo ""
echo "Checking for circular symlinks..."
BROKEN=$(find "$DST" -type l ! -exec test -e {} \; -print 2>/dev/null || true)
if [ -z "$BROKEN" ]; then
  echo "  PASS: No broken/circular symlinks"
else
  echo "  FAIL: Found broken symlinks:"
  echo "$BROKEN"
fi

echo ""
echo "Checking excluded directories..."
for d in ".claude" "docs"; do
  if [ -d "$DST/$d" ]; then
    echo "  FAIL: $d/ still present"
  else
    echo "  PASS: $d/ excluded"
  fi
done

echo ""
echo "Done. The public release is at: $DST"
echo "You can now: cd $DST && git init && git add -A && git commit"
