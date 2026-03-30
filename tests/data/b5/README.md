# B5 Benchmark Reference Data

**Status**: COMPLETE - Data generated 2026-02-04

## System Definition

- **QM Region**: 1 water molecule (3 atoms: O, H, H) - Residue 395
- **MM Region**: 883 TIP3P water molecules (2649 atoms)
- **Total**: 884 water molecules (2652 atoms)
- **Box**: 3.0 x 3.0 x 3.0 nm cubic
- **Purpose**: QM/MM stability and PBC validation (Milestone M2)

## Generated Files

| File | Status | Description |
|---|---|---|
| `b5_initial.gro` | COMPLETE | Equilibrated structure (NVT 100ps at 300K) |
| `b5_initial.ndx` | COMPLETE | Index file with QM_water group |
| `b5_topol.top` | COMPLETE | GROMACS topology (AMBER99SB-ILDN/TIP3P) |
| `qm_water_isolated.gen` | COMPLETE | QM water geometry in DFTB+ gen format |
| `dftb_in.hsd` | COMPLETE | DFTB+ input (DFTB2/mio-1-1) |
| `dftb_template.hsd` | DRAFT | DFTB+ template for DFTB3/3ob-3-1 (SK files needed) |
| `provenance.json` | COMPLETE | Complete provenance record |

## QM Region Selection

The QM water molecule was selected as the one closest to the box center:
- **Residue number**: 395
- **Atom indices** (1-based): 1183 (O), 1184 (H), 1185 (H)
- **Distance to center**: 0.225 nm
- **Oxygen position**: (1.435, 1.677, 1.623) nm

## Reference DFTB+ Calculation

Method: DFTB2 (SCC-DFTB) with mio-1-1 parameters

| Property | Value | Unit |
|---|---|---|
| Total energy | -3.7957527543 | Hartree |
| Mermin free energy | -3.7957527543 | Hartree |
| Mulliken charge (O) | -0.317 | e |
| Mulliken charge (H1) | +0.159 | e |
| Mulliken charge (H2) | +0.159 | e |
| Dipole moment | 1.695 | Debye |
| SCC iterations | 9 | - |

**Note**: For publication, regenerate with DFTB3/3ob-3-1 when SK files are available.

## Generation Protocol

Followed: `docs/theory/US-041/b5_generation_protocol.md`

### Steps Performed:
1. Created water box with `gmx solvate` (884 waters)
2. Energy minimization (converged in 165 steps, Fmax = 90.06 kJ/mol/nm)
3. NVT equilibration (100 ps at 300 K, seed=12345)
4. Selected QM water closest to box center
5. Ran standalone DFTB+ reference calculation

## GOLDEN RULE COMPLIANCE

Per `CLAUDE.md`: All numerical values in this benchmark come from actual GROMACS and DFTB+ calculations performed on 2026-02-04. No values were fabricated.

Evidence files:
- `preparation/em.log` - Energy minimization output
- `preparation/nvt.log` - NVT equilibration output
- `qm_water_reference.log` - DFTB+ calculation output
- `detailed.out` - DFTB+ detailed output
- `results.tag` - DFTB+ machine-readable results

---

*Generated: 2026-02-04*
*Verified by: QA Validation Specialist (Agent 09)*
