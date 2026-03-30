# GRO-DFTB

An open-source QM/MM interface coupling [DFTB+](https://dftbplus.org/) with [GROMACS](https://www.gromacs.org/) for energy-conserving ground-state molecular dynamics.

GRO-DFTB integrates the self-consistent-charge density-functional tight-binding (SCC-DFTB) method into GROMACS through its native `IForceProvider`/`MDModule` architecture, enabling nanosecond-scale QM/MM simulations of solvated systems under periodic boundary conditions.

**Paper**: V. Vchirawongkwin, "GRO-DFTB: Integrating SCC-DFTB with GROMACS for Energy-Conserving QM/MM Molecular Dynamics" (*submitted*)

## Features

- **C-ABI driver library** (`libgrodftb.so`) wrapping the DFTB+ shared-library API for in-process QM evaluation without file I/O
- **Two electrostatic embedding modes**:
  - Gaussian-damped cutoff with a C²-continuous quintic switching function
  - PME-compatible embedding that inserts self-consistent Mulliken charges into the GROMACS Ewald solver
- **Link-atom force projection** for covalent QM/MM boundaries with energy-consistent chain-rule redistribution
- **Bitwise-reproducible checkpoint/restart** of the QM charge state via GROMACS key-value-tree serialization
- **NVE energy conservation** validated on six benchmark systems (3–32 QM atoms), with normalized drift rates of −4.5 × 10⁻⁵ (cutoff) and (6.1 ± 1.0) × 10⁻⁴ kJ mol⁻¹ ps⁻¹ atom⁻¹ (PME, 10 seeds), both below the 0.005 acceptance threshold
- Throughput: 0.8–11 ns/day (PME mode, 8 OpenMP threads, 3–32 QM atoms)

## Architecture

```
Layer 5:  User Workflow (GROMACS mdrun)
Layer 4:  GROMACS QM/MM Backend         [Module 3]
Layer 3:  Embedding Engine [Module 4]   Link Atoms [Module 5]
Layer 2:  Unit Conversion               [Module 2]
Layer 1:  DFTB+ Driver Library          [Module 1]  (libgrodftb.so)
Layer 0:  DFTB+ (libdftbplus.so)        GROMACS (libgromacs.so)
```

Modules 1–5 are implemented in this release. Modules 6–9 (charge-transfer CVs, coupled-perturbed response, excited-state dynamics) are planned for future releases.

## Requirements

| Dependency | Version | Notes |
|---|---|---|
| DFTB+ | 25.1 (commit `fd31d873`) | Build with `WITH_API=ON`, `BUILD_SHARED_LIBS=ON` |
| GROMACS | 2027.0-dev (tag `v2026.0`) | Build with `GMX_DFTB=ON` |
| CMake | ≥ 3.28 | |
| C compiler | C11 (GCC ≥ 11 or Clang ≥ 14) | |
| C++ compiler | C++17 | Required for GROMACS backend |
| Fortran compiler | Any | Required for linking DFTB+ |

## Building

```bash
# 1. Clone with submodules
git clone --recurse-submodules https://github.com/ExaPsi/GRO-DFTB.git
cd GRO-DFTB

# 2. Build DFTB+
bash tools/build_dftbplus.sh

# 3. Build libgrodftb
cmake -B build \
  -DDFTBPlus_ROOT=external/dftbplus/_install
cmake --build build

# 4. Run tests
ctest --test-dir build

# 5. Build GROMACS with DFTB+ support
bash tools/build_gromacs.sh
```

See `examples/b5_tutorial/README.md` for a complete tutorial on running a QM/MM NVE trajectory.

## Benchmark Systems

| ID | System | N_QM | N_total | Purpose |
|---|---|---|---|---|
| B1 | Water dimer (gas) | 6 | 6 | Driver validation |
| B2 | Formamide (gas) | 6 | 6 | SCC convergence |
| B4 | Water dimer + 1 charge | 6 | 7 | Embedding validation |
| B5 | QM water in TIP3P box | 3 | 2652 | NVE stability (cutoff + PME) |
| B6 | Alanine dipeptide | 12 | 3163 | Link atoms + solvated NVE |
| — | Alanine tripeptide | 32 | 4271 | QM scaling |

All benchmark input files and reference data with full provenance are in `tests/data/`.

## Validation Summary

All analytic forces validated against central finite differences (δ = 10⁻⁴ Bohr, combined tolerance max(10⁻⁴|F|, 10⁻⁶ Ha/Bohr)).

| Test | System | Criterion | Achieved | Margin |
|---|---|---|---|---|
| Gas-phase energy match | B1, B2 | < 10⁻⁸ Ha | < 10⁻¹⁴ Ha | > 10⁶× |
| FD QM forces | B1 | < 10⁻⁴ rel | 2.81 × 10⁻⁶ | 36× |
| FD MM back-reaction | B4 | < 10⁻⁴ rel | 6.15 × 10⁻⁸ | 1600× |
| NVE drift (cutoff) | B5 | < 0.005 | −4.53 × 10⁻⁵ | 110× |
| NVE drift (PME) | B5 | < 0.005 | (6.1 ± 1.0) × 10⁻⁴ | 8.1× |
| NVE drift (link atoms) | B6 | < 0.005 | (6.3 ± 1.5) × 10⁻⁴ | 7.9× |
| NVE drift (32 QM) | Tripeptide | < 0.005 | (1.55 ± 0.35) × 10⁻³ | 3.2× |
| Checkpoint restart | B5 | |ΔE| < 10⁻⁶ Ha | 0.0 Ha | exact |

Drift units: kJ mol⁻¹ ps⁻¹ atom⁻¹. Multi-seed results report mean ± std (5–10 independent trajectories).

## Repository Structure

```
CMakeLists.txt          Top-level build
cmake/                  FindDFTBPlus.cmake, FindGROMACS.cmake
include/grodftb/        Public C headers (driver, embedding, switching, damping, linkatom, units, error)
src/                    Implementation (C11)
  driver/               DFTB+ driver library (grodftb_init, compute, finalize)
  embedding/            Switching function, Gaussian damping
  linkatom/             Link atom placement and force projection
  units/                CODATA 2022 unit conversion constants
tests/                  39 automated tests (30 libgrodftb + 9 GROMACS GTests)
  data/                 Benchmark inputs, SK parameters, provenance metadata
tools/                  Analysis scripts, build helpers, CI runner
examples/b5_tutorial/   Complete QM/MM NVE tutorial
```

## Analysis Tools

| Script | Purpose |
|---|---|
| `tools/analyze_pme_drift.py` | NVE drift analysis with block decomposition |
| `tools/analyze_nstlist_sweep.py` | Pair-search interval sensitivity |
| `tools/analyze_pbc_regression.py` | Box-size independence test |
| `tools/analyze_sigma_sweep.py` | Damping parameter sensitivity |
| `tools/analyze_energy_decomposition.py` | Energy component time series |
| `tools/analyze_a2_statistical.py` | Multi-seed statistical analysis |
| `tools/analyze_b6_multiseed.py` | Solvated B6 multi-seed analysis |
| `tools/parse_charges_bin.py` | DFTB+ charges.bin parser |
| `tools/compute_delta_pbc.py` | PBC correction residual analysis |

## Known Limitations

- **NVE only**: The QM-internal virial is not computed; NPT simulations are not validated.
- **DFTB2 only in dynamics**: DFTB3 energies are reproduced, but an upstream DFTB+ C API bug prevents gradient extraction with `ThirdOrderFull`.
- **Single-molecule QM regions**: Multi-molecule QM regions are limited by intermolecular double-counting in the charge-inserted Ewald approach.
- **Cutoff embedding**: Susceptible to a polarization instability in condensed-phase systems due to missing long-range screening; use PME embedding for production.
- **Static QM region**: Adaptive QM/MM boundary methods are not supported.

## Planned Extensions

1. **Coupled-perturbed SCC-DFTB** (dq/dR) for force-consistent biasing along charge-transfer collective variables
2. **TD-DFTB excited-state interfaces** for nonadiabatic dynamics via surface hopping or Ehrenfest propagation
3. **Intermolecular QM–QM exclusions** for stable multi-molecule QM regions

## Citation

If you use GRO-DFTB, please cite:

```bibtex
@article{Vchirawongkwin2026,
  author  = {Vchirawongkwin, Viwat},
  title   = {{GRO-DFTB}: Integrating {SCC-DFTB} with {GROMACS} for
             Energy-Conserving {QM/MM} Molecular Dynamics},
  year    = {2026},
  note    = {Manuscript submitted}
}
```

See also [CITATION.cff](CITATION.cff) for machine-readable citation metadata.

## License

[LGPL-3.0-or-later](LICENSE)

## Acknowledgments

The author thanks Associate Professor Dr. Chinapong Kritayakornupong (King Mongkut's University of Technology Thonburi) for discussions, the Department of Chemistry, Faculty of Science, Chulalongkorn University for research facilities, and the DFTB+ development team for maintaining the stable shared-library API.
