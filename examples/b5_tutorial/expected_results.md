# Expected Results for B5 QM/MM Tutorial

All values below come from **real calculations** archived in the GRO-DFTB repository.
No values are fabricated, estimated, or assumed.

## Source

**Story**: US-043b (Gaussian damping + PBC fix)
**Provenance file**: `tests/data/b5/drift_analysis_damped_500ps.json`
**Trajectory files**: `tests/data/b5/nve_damped/`

**Simulation parameters**:
- DFTB+ 25.1 (commit fd31d873) with mio-1-1 Slater-Koster parameters
- GROMACS 2027.0-dev (commit e1913f698e) with GMX_DFTB=ON
- Cutoff embedding: 1.2 nm with 0.2 nm quintic switching
- Gaussian damping: sigma = 0.1 nm
- SCCTolerance = 1e-10, Anderson mixer, Fermi filling 100 K
- Timestep: 0.5 fs, NVE ensemble, 500 ps total

## 500 ps Production Results

| Quantity | Value | Threshold | Status |
|---|---|---|---|
| NVE drift | -4.53 x 10^-5 kJ/mol/ps/atom | < 0.01 | PASS (221x margin) |
| Total energy mean | -38982.79 kJ/mol | — | — |
| Total energy sigma | 17.36 kJ/mol | < 450 | PASS |
| Temperature mean | 306.26 K | 300 +/- 10 | PASS |
| Temperature sigma | 4.77 K | < 10 | PASS |
| QM energy mean | -10704.68 kJ/mol | — | — |
| QM energy sigma | 0.92 kJ/mol | < 100 | PASS |
| Anomalous QM frames | 0 | 0 | PASS |
| Duration | 500.0 ps | — | — |
| Total frames | 100,001 | — | — |

## 10 ps Quick-Validation Results

For the 10 ps quick-validation run, expect:
- Drift: same order of magnitude (< 0.01 threshold easily satisfied)
- Temperature: within a few K of the 500 ps mean
- QM energy: consistent with the 500 ps range
- Zero SCC anomalies

Exact 10 ps values depend on starting conditions but should match the 500 ps
statistics within statistical noise.

## How These Values Were Obtained

1. B5 system equilibrated via 100 ps NVT at 300 K (V-rescale thermostat)
2. Production NVE run: 1,000,000 steps at 0.5 fs = 500 ps
3. Energy extracted with `gmx energy` (Total-Energy, Temperature, Quantum-En.)
4. Drift computed by linear regression of total energy vs. time
5. Results archived in `tests/data/b5/drift_analysis_damped_500ps.json`
6. GROMACS log confirms clean completion: `tests/data/b5/nve_damped/b5_nve_damped.log`
