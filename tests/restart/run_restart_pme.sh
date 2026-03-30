#!/usr/bin/env bash
# US-052 Restart Validation: PME Embedding Mode
#
# ThDD:Eq. T-US-052-N-1.1 -- trajectory comparison protocol
#
# Runs a B5 QM/MM simulation (1 QM water in 883 TIP3P waters) for 40 steps
# in two ways:
#   1. Continuous: 40 steps uninterrupted
#   2. Restarted: 20 steps -> checkpoint -> restart -> 20 more steps
#
# Compares QM energies at each post-restart step (21-40).
# Tolerance: 1e-6 Ha (ThDD:Eq. T-US-052-V-1.2)
#
# Prerequisites:
#   - GROMACS built with GMX_DFTB=ON at external/gromacs/_install_dftb/
#   - B5 system data at tests/data/b5/
#   - SK parameters at tests/data/slako/mio-1-1/
#
# Usage:
#   bash tests/restart/run_restart_pme.sh
#   GMX=/path/to/gmx bash tests/restart/run_restart_pme.sh
#
# LGPL-3.0-or-later

set -euo pipefail

# --- Configuration ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

GMX="${GMX:-${REPO_ROOT}/external/gromacs/_install_dftb/bin/gmx}"
NSTEPS_TOTAL=40
NSTEPS_SEG1=20
RESTART_STEP=21  # First post-restart step

# Source files
MDP_SRC="${REPO_ROOT}/examples/b5_tutorial/nve_qmmm_pme.mdp"
HSD_SRC="${REPO_ROOT}/examples/b5_tutorial/dftb_in.hsd"
GEN_SRC="${REPO_ROOT}/tests/data/b5/qm_water_isolated.gen"
GRO_SRC="${REPO_ROOT}/tests/data/b5/b5_initial.gro"
TOP_SRC="${REPO_ROOT}/tests/data/b5/b5_topol.top"
NDX_SRC="${REPO_ROOT}/tests/data/b5/b5_initial.ndx"
SLAKO_SRC="${REPO_ROOT}/tests/data/slako"

COMPARE_SCRIPT="${SCRIPT_DIR}/compare_energies.py"
ARCHIVE_DIR="${REPO_ROOT}/tests/data/b5_restart/pme"

echo "========================================================================"
echo "US-052 Restart Validation: PME Embedding Mode"
echo "========================================================================"
echo "GMX binary: ${GMX}"
echo "Repository: ${REPO_ROOT}"
echo "Total steps: ${NSTEPS_TOTAL}"
echo "Segment 1 steps: ${NSTEPS_SEG1}"
echo ""

# --- Validate prerequisites ---
for f in "${GMX}" "${MDP_SRC}" "${HSD_SRC}" "${GEN_SRC}" "${GRO_SRC}" "${TOP_SRC}" "${NDX_SRC}" "${COMPARE_SCRIPT}"; do
    if [ ! -f "${f}" ]; then
        echo "ERROR: Required file not found: ${f}" >&2
        exit 1
    fi
done
if [ ! -d "${SLAKO_SRC}" ]; then
    echo "ERROR: SK parameter directory not found: ${SLAKO_SRC}" >&2
    exit 1
fi

# --- Create temporary working directory ---
WORKDIR="$(mktemp -d /tmp/us052_pme_XXXXXX)"
echo "Working directory: ${WORKDIR}"

cleanup() {
    echo ""
    echo "Temporary files in: ${WORKDIR}"
    echo "(Remove manually with: rm -rf ${WORKDIR})"
}
trap cleanup EXIT

# --- Create test MDP with nsteps=40, nstenergy=1 ---
TEST_MDP="${WORKDIR}/test.mdp"
cp "${MDP_SRC}" "${TEST_MDP}"
# Override nsteps to 40 and nstenergy to 1 (every step)
sed -i 's/^nsteps\b.*/nsteps                   = '"${NSTEPS_TOTAL}"'/' "${TEST_MDP}"
sed -i 's/^nstenergy\b.*/nstenergy                = 1/' "${TEST_MDP}"

echo "Test MDP created with nsteps=${NSTEPS_TOTAL}, nstenergy=1"

# --- Run grompp (once, shared tpr) ---
echo ""
echo "--- Running grompp ---"
TPR="${WORKDIR}/topol.tpr"
${GMX} grompp \
    -f "${TEST_MDP}" \
    -c "${GRO_SRC}" \
    -p "${TOP_SRC}" \
    -n "${NDX_SRC}" \
    -o "${TPR}" \
    -maxwarn 1 \
    2>&1 | tail -5

if [ ! -f "${TPR}" ]; then
    echo "ERROR: grompp failed to produce tpr file" >&2
    exit 1
fi
echo "tpr created: ${TPR}"

# --- Set up continuous run directory ---
CONT_DIR="${WORKDIR}/cont"
mkdir -p "${CONT_DIR}"
cp "${TPR}" "${CONT_DIR}/topol.tpr"
cp "${HSD_SRC}" "${CONT_DIR}/dftb_in.hsd"
cp "${GEN_SRC}" "${CONT_DIR}/qm_water_isolated.gen"
ln -s "${SLAKO_SRC}" "${CONT_DIR}/slako"

# --- Set up restart run directory ---
RESTART_DIR="${WORKDIR}/restart"
mkdir -p "${RESTART_DIR}"
cp "${TPR}" "${RESTART_DIR}/topol.tpr"
cp "${HSD_SRC}" "${RESTART_DIR}/dftb_in.hsd"
cp "${GEN_SRC}" "${RESTART_DIR}/qm_water_isolated.gen"
ln -s "${SLAKO_SRC}" "${RESTART_DIR}/slako"

# --- Continuous run: 40 steps ---
echo ""
echo "--- Continuous run (${NSTEPS_TOTAL} steps) ---"
cd "${CONT_DIR}"
${GMX} mdrun \
    -s topol.tpr \
    -deffnm cont \
    -pin off \
    -ntmpi 1 \
    -ntomp 1 \
    2>&1 | tail -3
echo "Continuous run complete."

if [ ! -f "${CONT_DIR}/cont.edr" ]; then
    echo "ERROR: Continuous run failed — no .edr file" >&2
    exit 1
fi

# --- Restart run: Segment 1 (20 steps) ---
echo ""
echo "--- Restart run: Segment 1 (${NSTEPS_SEG1} steps) ---"
cd "${RESTART_DIR}"
${GMX} mdrun \
    -s topol.tpr \
    -deffnm restart \
    -pin off \
    -ntmpi 1 \
    -ntomp 1 \
    -nsteps ${NSTEPS_SEG1} \
    2>&1 | tail -3
echo "Segment 1 complete."

if [ ! -f "${RESTART_DIR}/restart.cpt" ]; then
    echo "ERROR: Segment 1 failed — no checkpoint file" >&2
    exit 1
fi

# --- Restart run: Segment 2 (continue from checkpoint) ---
echo ""
echo "--- Restart run: Segment 2 (continue from step ${NSTEPS_SEG1} to ${NSTEPS_TOTAL}) ---"
cd "${RESTART_DIR}"
${GMX} mdrun \
    -s topol.tpr \
    -deffnm restart \
    -pin off \
    -ntmpi 1 \
    -ntomp 1 \
    -cpi restart.cpt \
    2>&1 | tail -3
echo "Segment 2 complete."

# --- Check for SCC failures (AC-4) ---
echo ""
echo "--- Checking for SCC failures ---"
SCC_FAILURES=$(grep -c "SCC not converged" "${RESTART_DIR}/restart.log" 2>/dev/null || true)
if [ "${SCC_FAILURES}" -gt 0 ]; then
    echo "WARNING: ${SCC_FAILURES} SCC failure(s) detected in restart run"
else
    echo "No SCC failures detected (PASS)"
fi

# --- Compare energies ---
echo ""
echo "--- Comparing QM energies ---"
RESULTS_JSON="${WORKDIR}/results_pme.json"
cd "${REPO_ROOT}"
python3 "${COMPARE_SCRIPT}" \
    --gmx "${GMX}" \
    --cont-edr "${CONT_DIR}/cont.edr" \
    --restart-edr "${RESTART_DIR}/restart.edr" \
    --restart-step ${RESTART_STEP} \
    --output "${RESULTS_JSON}" \
    --mode pme

COMPARE_EXIT=$?

# --- Archive results ---
echo ""
echo "--- Archiving results ---"
mkdir -p "${ARCHIVE_DIR}"
cp "${RESULTS_JSON}" "${ARCHIVE_DIR}/results.json"
echo "Results archived to: ${ARCHIVE_DIR}/results.json"

# --- Final summary ---
echo ""
echo "========================================================================"
if [ ${COMPARE_EXIT} -eq 0 ] && [ "${SCC_FAILURES}" -eq 0 ]; then
    echo "US-052 PME RESTART VALIDATION: PASS"
    echo "  - QM energy match: PASS (within 1e-6 Ha)"
    echo "  - SCC convergence: PASS (0 failures)"
else
    echo "US-052 PME RESTART VALIDATION: FAIL"
    if [ ${COMPARE_EXIT} -ne 0 ]; then
        echo "  - QM energy match: FAIL"
    fi
    if [ "${SCC_FAILURES}" -gt 0 ]; then
        echo "  - SCC convergence: FAIL (${SCC_FAILURES} failures)"
    fi
fi
echo "========================================================================"

exit ${COMPARE_EXIT}
