# US-052 Restart Validation Tests

## Overview

These tests validate that the GRO-DFTB checkpoint/restart implementation
(US-051) produces bit-identical (or near-identical within SCC convergence
bounds) QM energy results when a B5 QM/MM simulation is interrupted and
restarted from a GROMACS checkpoint.

## Test Protocol

The protocol follows ThDD:Eq. T-US-052-N-1.1 (trajectory comparison):

1. **Continuous run**: Run B5 system for 40 steps without interruption
2. **Restarted run**: Run 20 steps, checkpoint, restart, run 20 more steps
3. **Compare**: Extract `Quantum-En.` from both `.edr` files and compute
   the energy difference at each post-restart step (21-40)

Both runs use the same `.tpr` file (identical initial conditions and
parameters). Energy is recorded at every step (`nstenergy = 1`).

## Acceptance Criteria

| AC | Criterion | Tolerance | Theory Ref |
|----|-----------|-----------|------------|
| AC-1/2 | QM energy match at step 21 | 1e-6 Ha | Eq. T-US-052-V-1.1/V-1.2 |
| AC-3 | Multi-step max |dE| (steps 21-40) | 1e-6 Ha | Eq. T-US-052-N-2.1 |
| AC-4 | No SCC failure on restart | 0 failures | Eq. T-US-052-V-2.1 |
| AC-5 | No energy discontinuity | within 3-sigma | Eq. T-US-052-V-2.2 |

Expected result: bitwise identity (|dE| = 0.0 Ha) based on US-051
empirical evidence.

## Prerequisites

- GROMACS built with `GMX_DFTB=ON` at `external/gromacs/_install_dftb/bin/gmx`
- B5 system data at `tests/data/b5/`
- Slater-Koster parameters at `tests/data/slako/mio-1-1/`
- Python 3 (for `compare_energies.py`)

## Files

| File | Description |
|------|-------------|
| `run_restart_cutoff.sh` | Cutoff embedding restart test |
| `run_restart_pme.sh` | PME embedding restart test |
| `compare_energies.py` | Energy comparison and analysis script |
| `README.md` | This documentation |

## How to Run

### Cutoff mode
```bash
bash tests/restart/run_restart_cutoff.sh
```

### PME mode
```bash
bash tests/restart/run_restart_pme.sh
```

### Custom GMX binary
```bash
GMX=/path/to/custom/gmx bash tests/restart/run_restart_cutoff.sh
```

## Output

Each test writes results to:
- `tests/data/b5_restart/cutoff/results.json` (cutoff mode)
- `tests/data/b5_restart/pme/results.json` (PME mode)

The JSON file contains:
- Per-step energy differences in Hartree and kJ/mol
- Maximum absolute deviation and RMS deviation
- Energy discontinuity analysis at the restart point
- Full provenance (file paths, theory references)

## Theory References

| Equation | Description |
|----------|-------------|
| Eq. T-US-052-N-1.1 | Trajectory comparison protocol |
| Eq. T-US-052-N-2.1 | Maximum absolute energy deviation |
| Eq. T-US-052-N-2.2 | RMS energy deviation |
| Eq. T-US-052-U-1.2 | kJ/mol to Hartree conversion (1 Ha = 2625.5 kJ/mol) |
| Eq. T-US-052-V-1.1 | Cutoff mode tolerance: 1e-6 Ha |
| Eq. T-US-052-V-1.2 | PME mode tolerance: 1e-6 Ha |
| Eq. T-US-052-V-2.1 | SCC convergence on restart |
| Eq. T-US-052-V-2.2 | Energy discontinuity within 3-sigma |

## Determinism Notes

- Both runs use `-ntmpi 1 -ntomp 1` for deterministic comparison (R-2 mitigation)
- `-pin off` is always used (avoids thread affinity issues in scripted execution)
- DFTB+ resolves paths from CWD; each run directory contains dftb_in.hsd,
  qm_water_isolated.gen, and a symlink to slako/
