# B4 Damped Reference Data Provenance

## Calculation Details
- **System**: Water dimer (6 atoms) + 1 MM point charge at (10, 0, 0) Bohr
- **Method**: SCC-DFTB with mio-1-1 Slater-Koster parameters
- **Embedding**: GaussianBlurWidth [Bohr] = 3.78
- **MM charge**: +0.417 e
- **Software**: DFTB+ 25.1 (commit fd31d873)
- **Date**: 2026-02-08

## Input
- `dftb_in.hsd` — standalone DFTB+ input with PointCharges and GaussianBlurWidth
- Geometry from `tests/data/b4/geo.gen` (water dimer in Bohr)

## Reference Values (from autotest.tag)
- **mermin_energy**: -8.15593240308962 Ha
- **forces**: 18 components (6 atoms x 3)
- **forces_ext_charges**: 3 components (1 charge x 3)

## Reference Values (from detailed.out)
- **Gross charges**: [-0.59496640, 0.29861031, 0.29861031, -0.67376910, 0.34884831, 0.32266657]
- **SCC converged**: 9 iterations, error = 6.99e-11

## Used In
- `tests/test_damping_integration.c`: AC-10 (test_b4_damped_energy_match)
- Energy tolerance: 1e-8 Ha
- Force tolerance: 1e-8 Ha/Bohr
- Charge tolerance: 1e-6 e

## Notes
- This directory is a protected copy of the standalone DFTB+ calculation outputs.
- The parent directory (`tests/data/b4_damped/`) contains the API HSD (`dftb_in.hsd`)
  used by the test, which writes its own output files there (overwriting any existing ones).
- Reference outputs must be kept in this `reference/` subdirectory to avoid being overwritten.
