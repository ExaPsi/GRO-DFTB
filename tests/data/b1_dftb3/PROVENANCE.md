# B1 DFTB3/3OB Water Dimer Reference Data

## System
- Water dimer (H2O)2, gas phase, 6 atoms
- Geometry from tests/data/b1/geo.gen (same as DFTB2 validation)

## Method
- DFTB3 (ThirdOrderFull = Yes)
- 3OB-3-1 Slater-Koster parameters
- HubbardDerivs: O = -0.1575, H = -0.1857
- SCCTolerance = 1.0e-10
- Broyden mixer

## Software
- DFTB+ 25.1 (commit fd31d873)
- Built with WITH_API=ON, BUILD_SHARED_LIBS=ON

## Results (standalone dftb+)
- Mermin free energy: -8.1286981909 Ha (-221.1931 eV)
- SCC iterations: 11
- Forces: see results.tag
- Mulliken charges: O1=-0.6738, H2=+0.3204, H3=+0.3366,
                     O4=-0.6399, H5=+0.3284, H6=+0.3284

## API Validation
- Energy via C API: -8.1286981909 Ha (error: 3.79e-12 Ha)
- Charges via C API: match standalone
- **Gradients: NOT AVAILABLE** — dftbp_get_gradients() causes Fortran
  ERROR STOP with ThirdOrderFull in DFTB+ 25.1 C API
- This is a DFTB+ API limitation, not a GRO-DFTB bug

## DFTB2 vs DFTB3 Comparison
- DFTB2/mio-1-1 energy: -4.0760783640 Ha (from tests/data/b1/)
- DFTB3/3ob-3-1 energy: -8.1286981909 Ha
- Difference reflects different parametrization, not an error

## Date
- 2026-02-12
