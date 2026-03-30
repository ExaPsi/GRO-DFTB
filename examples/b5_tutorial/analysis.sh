#!/bin/bash
# GRO-DFTB Tutorial: B5 NVE Analysis Script
#
# Usage: ./analysis.sh [EDR_FILE] [N_ATOMS]
#
# Extracts energy, temperature, and QM energy from a GROMACS .edr file,
# computes NVE drift, and prints a pass/fail validation checklist.
#
# Requirements: gmx (GROMACS), python3 with numpy
#
# LGPL-3.0-or-later

set -euo pipefail

EDR="${1:-b5_qmmm.edr}"
NATOMS="${2:-2652}"
GMX="${GMX:-gmx}"

echo "============================================================"
echo "  GRO-DFTB B5 NVE Analysis"
echo "============================================================"
echo "  EDR file:  ${EDR}"
echo "  N_atoms:   ${NATOMS}"
echo ""

if [ ! -f "${EDR}" ]; then
    echo "ERROR: EDR file '${EDR}' not found."
    echo "Run gmx mdrun first, then re-run this script."
    exit 1
fi

# --- Step 1: Extract energy time series ---
echo ">>> Extracting Total-Energy..."
echo "Total-Energy" | ${GMX} energy -f "${EDR}" -o total_energy.xvg 2>/dev/null

echo ">>> Extracting Temperature..."
echo "Temperature" | ${GMX} energy -f "${EDR}" -o temperature.xvg 2>/dev/null

echo ">>> Extracting Quantum Energy..."
echo "Quantum-En." | ${GMX} energy -f "${EDR}" -o quantum_energy.xvg 2>/dev/null

echo ""

# --- Step 2: Compute NVE drift and statistics ---
echo ">>> Computing drift and statistics..."
python3 - "${NATOMS}" <<'PYEOF'
import sys
import numpy as np

natoms = int(sys.argv[1])

def load_xvg(filename):
    """Load a GROMACS .xvg file, skipping comment/header lines."""
    t, v = [], []
    with open(filename) as f:
        for line in f:
            if line.startswith(('#', '@')):
                continue
            parts = line.split()
            if len(parts) >= 2:
                t.append(float(parts[0]))
                v.append(float(parts[1]))
    return np.array(t), np.array(v)

# Total energy
t_e, etot = load_xvg("total_energy.xvg")
slope, intercept = np.polyfit(t_e, etot, 1)
drift = abs(slope) / natoms
e_mean = np.mean(etot)
e_std = np.std(etot)

# Temperature
t_t, temp = load_xvg("temperature.xvg")
t_mean = np.mean(temp)
t_std = np.std(temp)

# QM energy
t_q, eqm = load_xvg("quantum_energy.xvg")
qm_mean = np.mean(eqm)
qm_std = np.std(eqm)
qm_threshold = 5.0  # kJ/mol for anomaly detection
anomalous = np.sum(np.abs(eqm - qm_mean) > qm_threshold)

duration = t_e[-1] - t_e[0]

print("")
print("============================================================")
print("  Results")
print("============================================================")
print(f"  Duration:             {duration:.1f} ps ({len(etot)} frames)")
print(f"  Total energy:         {e_mean:.2f} +/- {e_std:.2f} kJ/mol")
print(f"  NVE drift:            {drift:.2e} kJ/mol/ps/atom")
print(f"  Temperature:          {t_mean:.2f} +/- {t_std:.2f} K")
print(f"  QM energy:            {qm_mean:.2f} +/- {qm_std:.2f} kJ/mol")
print(f"  Anomalous QM frames:  {int(anomalous)}")
print("")

# --- Validation checklist ---
DRIFT_THRESHOLD = 0.01       # kJ/mol/ps/atom (specs.md S21.1)
TEMP_MEAN_LOW = 290.0        # K
TEMP_MEAN_HIGH = 310.0       # K
TEMP_STD_THRESHOLD = 10.0    # K
QM_STD_THRESHOLD = 100.0     # kJ/mol

checks = [
    ("NVE drift < 0.01 kJ/mol/ps/atom",  drift < DRIFT_THRESHOLD),
    ("Temperature 290-310 K",             TEMP_MEAN_LOW < t_mean < TEMP_MEAN_HIGH),
    ("Temperature sigma < 10 K",          t_std < TEMP_STD_THRESHOLD),
    ("QM energy sigma < 100 kJ/mol",      qm_std < QM_STD_THRESHOLD),
    ("Zero SCC anomalies",                anomalous == 0),
]

print("============================================================")
print("  Validation Checklist")
print("============================================================")
all_pass = True
for label, passed in checks:
    status = "PASS" if passed else "FAIL"
    marker = "[x]" if passed else "[ ]"
    if not passed:
        all_pass = False
    print(f"  {marker} {label}: {status}")

print("")
if all_pass:
    print("  OVERALL: PASS -- All checks passed.")
else:
    print("  OVERALL: FAIL -- One or more checks failed. See above.")
print("============================================================")
PYEOF
