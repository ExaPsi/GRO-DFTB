# Changelog

All notable changes to GRO-DFTB are documented in this file.

Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html)

## [v1.0] — 2026

Paper release accompanying: V. Vchirawongkwin, "GRO-DFTB: Integrating SCC-DFTB with GROMACS for Energy-Conserving QM/MM Molecular Dynamics" (submitted).

### Highlights
- Complete ground-state QM/MM infrastructure for SCC-DFTB in GROMACS
- Two electrostatic embedding modes: Gaussian-damped cutoff and PME-compatible
- Link-atom force projection for covalent QM/MM boundaries
- NVE energy conservation validated on six benchmark systems (3–32 QM atoms)
- Bitwise-reproducible checkpoint/restart
- Throughput: 0.8–11 ns/day (PME, 8 threads)

## [v0.2] — 2026-02-10

M2b: PME-compatible embedding.

### Added
- PME-compatible electrostatic embedding via charge-inserted Ewald approach
- Pair-search-synchronized Mulliken charge insertion into GROMACS PME solver
- Real-space erfc(αr)/r kernel with Gaussian damping combination
- Reciprocal-space potential and field extraction via B-spline interpolation
- Energy double-counting correction for QM charges in the Ewald sum
- Runtime embedding consistency checks (three algebraic relations)
- Checkpoint/restart of QM charge state via GROMACS KeyValueTree serialization
- PBC regression tests across three box sizes (2.5, 3.0, 4.0 nm)
- Damping parameter sensitivity sweep (σ = 0.05–0.20 nm)
- 10-seed PME NVE drift analysis with statistical reproducibility

### Validated
- PME NVE drift: (6.1 ± 1.0) × 10⁻⁴ kJ/mol/ps/atom (10 seeds, 8.1× margin)
- Solvated B6 with link atoms: (6.3 ± 1.5) × 10⁻⁴ (5 seeds)
- Checkpoint restart: |ΔE| = 0.0 Ha (bitwise identical, both modes)
- PME erfc kernel FD force error: 3.66 × 10⁻¹⁰

## [v0.1] — 2026-02-02

M1: DFTB+ driver library + cutoff QM/MM coupling.

### Added
- `libgrodftb.so`: C-ABI driver library wrapping DFTB+ shared-library API
- Handle-based lifecycle: init → set geometry → set embedding → compute → retrieve → finalize
- Gradient negation (DFTB+ returns +dE/dR; library returns forces −dE/dR)
- CODATA 2022 unit conversion (Bohr↔nm, Hartree↔kJ/mol) with dual-constant policy
- `DFTBForceProvider` implementing GROMACS `IForceProvider` interface
- Cutoff electrostatic embedding with Gaussian damping (erf(r/σ)/r)
- C²-continuous quintic switching function for smooth cutoff truncation
- Link-atom placement and energy-consistent force projection
- Charge redistribution schemes (none/zero/shift) for QM/MM boundary
- PBC-aware coordinate unwrapping for contiguous QM clusters
- GROMACS QM/MM preprocessing reuse (topology splitting, charge zeroing, NB exclusions)
- `GMX_DFTB` CMake flag with stub/implementation pattern
- 30 automated tests for libgrodftb + 9 GROMACS-side GTests
- Benchmark data for B1 (water dimer), B2 (formamide), B4 (embedded dimer), B5 (solvated water), B6 (alanine dipeptide)
- B5 tutorial with step-by-step instructions

### Validated
- Gas-phase energy match: < 10⁻¹⁴ Ha vs standalone DFTB+ (B1, B2)
- FD force consistency: all components < 10⁻⁴ relative error
- Cutoff NVE drift: −4.53 × 10⁻⁵ kJ/mol/ps/atom (500 ps, 110× margin)
- B6 link-atom force and torque conservation to machine precision
- Zero memory leaks (Valgrind + ASan/UBSan)
