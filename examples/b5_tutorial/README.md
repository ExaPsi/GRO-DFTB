# GRO-DFTB Tutorial: QM/MM Molecular Dynamics with DFTB+

This tutorial walks through a complete QM/MM molecular dynamics simulation
using GRO-DFTB — from input preparation to trajectory validation. By the end,
you will have run an energy-conserving NVE trajectory of a quantum-mechanical
water molecule embedded in a classical TIP3P solvent box.

**System**: B5 benchmark — 1 QM water (3 atoms, SCC-DFTB2) in 883 TIP3P
waters (2649 MM atoms), 3.0 nm cubic box.

**What you will learn**:
- How to set up DFTB+ input (HSD) and GROMACS parameters (MDP) for QM/MM
- How GROMACS preprocesses the topology for QM/MM
- How to run and monitor a QM/MM simulation
- How to validate energy conservation and detect common problems

---

## 1. Prerequisites

### Software

| Component | Version | Notes |
|---|---|---|
| GROMACS | 2027.0-dev or later | Built with `GMX_DFTB=ON` |
| DFTB+ | 25.1 | Built with `WITH_API=ON`, `BUILD_SHARED_LIBS=ON` |
| libgrodftb | 0.1.x | Linked into the GROMACS build |
| Python 3 | 3.8+ | With numpy (for analysis script) |

### Files

You need these files from the GRO-DFTB repository:

| File | Location | Description |
|---|---|---|
| Starting structure | `tests/data/b5/b5_initial.gro` | Pre-equilibrated coordinates with velocities |
| Index file | `tests/data/b5/b5_initial.ndx` | Defines `QM_water` atom group |
| Topology | `tests/data/b5/b5_topol.top` | AMBER99SB-ILDN + TIP3P |
| QM geometry | `tests/data/b5/qm_water_isolated.gen` | Initial QM atom positions (GenFormat) |
| SK parameters | `tests/data/slako/mio-1-1/` | DFTB2 Slater-Koster files for O and H |
| Force field | `external/gromacs/_install_dftb/share/gromacs/top/` | GROMACS force field directory |

### Working directory setup

Create a clean working directory and copy the required files:

```bash
# Create working directory
mkdir -p b5_run && cd b5_run

# Copy structural files from repository
cp /path/to/gro-dftb/tests/data/b5/b5_initial.gro .
cp /path/to/gro-dftb/tests/data/b5/b5_initial.ndx .
cp /path/to/gro-dftb/tests/data/b5/b5_topol.top .
cp /path/to/gro-dftb/tests/data/b5/qm_water_isolated.gen .

# Copy tutorial input files
cp /path/to/gro-dftb/examples/b5_tutorial/dftb_in.hsd .
cp /path/to/gro-dftb/examples/b5_tutorial/nve_qmmm.mdp .
cp /path/to/gro-dftb/examples/b5_tutorial/analysis.sh .

# Set up Slater-Koster parameters (choose one):
# Option A: Symlink (recommended)
mkdir -p slako
ln -s /path/to/gro-dftb/tests/data/slako/mio-1-1 slako/mio-1-1

# Option B: Copy
# mkdir -p slako/mio-1-1
# cp /path/to/gro-dftb/tests/data/slako/mio-1-1/*.skf slako/mio-1-1/

# Option C: Environment variable (advanced)
# export DFTBPLUS_PARAM_DIR=/path/to/slako
# Then edit dftb_in.hsd: change Prefix to "mio-1-1/"
```

> **Path resolution**: DFTB+ resolves all relative paths in the HSD file
> from the **current working directory**, not from the HSD file location.
> Always run `gmx mdrun` from the directory containing `dftb_in.hsd` and
> the `slako/` subdirectory.

---

## 2. System Overview

The B5 benchmark system consists of:

- **QM region**: 1 water molecule (3 atoms: O, H, H) treated with SCC-DFTB2
  using the mio-1-1 Slater-Koster parameter set
- **MM region**: 883 TIP3P water molecules (2649 atoms) treated with the
  AMBER99SB-ILDN force field
- **Box**: 3.0 nm cubic, periodic boundary conditions
- **Interactions**:
  - QM internal: DFTB+ (band energy + SCC Coulomb + repulsive potential)
  - QM-MM electrostatic: Cutoff embedding with quintic switching and Gaussian damping
  - QM-MM van der Waals: Classical Lennard-Jones (handled by GROMACS)
  - MM-MM: Full AMBER99SB-ILDN + PME electrostatics

There are no covalent bonds crossing the QM/MM boundary (the QM water is an
intact molecule), so no link atoms are needed.

---

## 3. Input Files

### 3.1 MDP Parameters

The MDP file (`nve_qmmm.mdp`) controls the GROMACS simulation. Key QM/MM
parameters:

| MDP Parameter | Value | Physical Justification |
|---|---|---|
| `integrator = md` | Leapfrog | Symplectic, time-reversible integrator |
| `dt = 0.0005` | 0.5 fs | ~20 samples per O-H stretch period (9.8 fs) |
| `nsteps = 20000` | 10 ps | Quick validation (set to 1000000 for 500 ps production) |
| `constraints = h-bonds` | LINCS | Constrains MM water O-H bonds only (see Section 5) |
| `tcoupl = no` | NVE | Microcanonical ensemble for drift testing |
| `nstenergy = 10` | Every 5 fs | Sufficient data points for drift analysis |
| `qmmm-cp2k-dftbplus-embedding = CUTOFF` | Direct sum | Correct for validation |
| `qmmm-cp2k-dftbplus-cutoff = 1.2` | nm | Outer cutoff for embedding charges |
| `qmmm-cp2k-dftbplus-switch-width = 0.2` | nm | Quintic switching region [1.0, 1.2] nm |
| `qmmm-cp2k-dftbplus-damping-sigma = 0.1` | nm | Gaussian damping of short-range Coulomb |
| `qmmm-cp2k-dftbplus-charge-redistribution = NONE` | — | No covalent QM/MM bonds in B5 |

**Why the `qmmm-cp2k-` prefix?** The DFTB+ backend is integrated into the
GROMACS QM/MM infrastructure that was originally built for CP2K. The
`qmmm-cp2k-` prefix is a naming convention inherited from this infrastructure,
not a dependency on CP2K software.

**Why `switch-width > 0`?** Without switching, the embedding interaction has a
hard cutoff — forces are discontinuous at 1.2 nm, causing energy drift of
1-10 kJ/mol/ps/atom. The quintic switching function S(u) = 1 - 10u^3 + 15u^4 - 6u^5
provides C2 continuity (continuous value, first, and second derivatives),
eliminating this artifact.

**Why `damping-sigma > 0`?** The bare Coulomb potential 1/r diverges when an
MM atom approaches a QM atom at short range. This can cause SCC convergence
failure and simulation crashes. Gaussian damping replaces 1/r with
erf(r/sigma)/r, which is finite at r=0. At distances above ~3*sigma (0.3 nm),
the damped potential recovers the bare Coulomb to better than 0.01%.

### 3.2 HSD Settings

The HSD file (`dftb_in.hsd`) configures DFTB+. Key settings for QM/MM
robustness:

| HSD Setting | Value | Why |
|---|---|---|
| `SCC = Yes` | — | Required for DFTB2 charge self-consistency |
| `SCCTolerance = 1e-10` | e | Tight convergence minimizes NVE drift |
| `MaxSCCIterations = 500` | — | Extra headroom for perturbed embedded systems |
| `ReadInitialCharges = No` | — | **CRITICAL** — see warning below |
| `ConvergentSCCOnly = No` | — | **CRITICAL** — see warning below |
| `Mixer = Anderson` | — | More robust than Broyden for QM/MM |
| `MixingParameter = 0.05` | — | Conservative mixing prevents oscillation |
| `Filling = Fermi { Temperature [K] = 100 }` | — | Smooths electron occupations |
| `ForceEvaluation = "dynamics"` | — | Pulay correction for Fermi filling |
| `ParserVersion = 14` | — | Parser version for DFTB+ 25.1 |
| `WriteDetailedOut = No` | — | Avoid per-step file I/O in MD hot loop |

> **WARNING: ReadInitialCharges MUST be No**
>
> When using DFTB+ through the GRO-DFTB API, charge persistence between MD
> steps is automatic (the DftbPlus object retains Mulliken charges in memory).
> Setting `ReadInitialCharges = Yes` causes DFTB+ to look for a `charges.bin`
> file on disk. If this file does not exist, DFTB+ triggers a Fortran
> `ERROR STOP` — an **unrecoverable process abort** that crashes GROMACS
> with no useful error message.

> **WARNING: ConvergentSCCOnly MUST be No**
>
> With `Yes`, a single SCC non-convergence event triggers Fortran `ERROR STOP`,
> crashing the simulation irrecoverably. With `No`, GRO-DFTB receives the
> unconverged result and the MD integrator continues. Occasional non-convergence
> during transient geometries is rare but possible in QM/MM.

---

## 4. What grompp Does to Your Topology

When you run `gmx grompp` with QM/MM enabled, GROMACS automatically modifies
the topology to prevent double-counting between the QM and MM force
calculations. Understanding these changes helps diagnose problems.

**Automatic preprocessing steps**:

1. **Charge zeroing**: Classical partial charges on QM atoms are set to zero.
   This prevents GROMACS PME from computing QM-MM electrostatics (which DFTB+
   handles via embedding).

2. **Nonbonded exclusions**: QM-QM nonbonded interactions (Coulomb + LJ) are
   excluded from the classical force calculation. DFTB+ is the sole provider
   of QM-QM forces.

3. **Bond conversion**: QM-internal bonds are converted to `F_CONNBOND`
   (connectivity-only, zero force). This preserves the topology graph for
   constraint detection but ensures no classical bond forces compete with
   DFTB+ quantum forces.

4. **Angle/dihedral removal**: Angles with 2+ QM atoms and dihedrals with 3+
   QM atoms are removed from the classical force calculation.

5. **Constraint check**: If QM atoms have constrained bonds, grompp issues a
   warning. For B5, the QM water's O-H bonds are converted to connectivity
   bonds in step 3, so LINCS has no effect on them — QM dynamics are fully
   governed by DFTB+.

**What is NOT changed**:
- **LJ parameters on QM atoms are retained**. QM-MM van der Waals interactions
  are handled classically by GROMACS. This provides the short-range repulsion
  that keeps MM atoms from overlapping with QM atoms.

---

## 5. Running grompp

Generate the run input file (.tpr):

```bash
gmx grompp \
    -f nve_qmmm.mdp \
    -c b5_initial.gro \
    -p b5_topol.top \
    -n b5_initial.ndx \
    -o b5_qmmm.tpr \
    -maxwarn 1
```

**Expected output**: grompp prints information about the QM/MM setup
(number of QM atoms, embedding method) and two informational NOTEs about
energy calculation frequency and NVE temperature. These are NOTEs (not
warnings) and are expected. The `-maxwarn 1` flag is a safety net for
future GROMACS versions.

**If grompp fails**, check:
- Is the `QM_water` group defined in `b5_initial.ndx`?
- Does the topology include the correct force field?
- Are all referenced files in the current directory?

---

## 6. Quick Validation (10 ps)

Before committing to a 20-hour production run, verify the setup with a short
10 ps trajectory (~2 minutes wall time):

```bash
gmx mdrun -deffnm b5_qmmm -pin off -ntmpi 1 -ntomp 4
```

**Required flags**:

| Flag | Why |
|---|---|
| `-deffnm b5_qmmm` | Consistent output file naming |
| `-pin off` | **Required** — prevents thread pinning issues (see Pitfalls) |
| `-ntmpi 1` | **Required** — QM/MM must run on a single MPI rank |
| `-ntomp 4` | OpenMP threads (adjust for your system; 4-8 recommended) |

**What to check after the 10 ps run**:

1. **Did it complete?** Check the log file:
   ```bash
   tail -5 b5_qmmm.log
   ```
   You should see "Finished mdrun" with the step count.

2. **Run the analysis script**:
   ```bash
   bash analysis.sh b5_qmmm.edr 2652
   ```
   All validation checks should show PASS.

3. **Is the total energy reasonable?** Should be around -39,000 kJ/mol for B5.

4. **Is the temperature reasonable?** Should be near 300 K (the NVT
   equilibration target).

If everything passes, proceed to the production run.

---

## 7. Production Run (500 ps)

For publication-quality drift analysis, run a longer trajectory. Edit
`nve_qmmm.mdp` and change:

```
nsteps = 1000000    ; 500 ps total
```

Then re-run grompp and launch the production run in the background:

```bash
# Re-generate .tpr with new nsteps
gmx grompp \
    -f nve_qmmm.mdp \
    -c b5_initial.gro \
    -p b5_topol.top \
    -n b5_initial.ndx \
    -o b5_qmmm.tpr \
    -maxwarn 1

# Launch in background (~20 hours wall time)
nohup gmx mdrun -deffnm b5_qmmm -pin off -ntmpi 1 -ntomp 8 &> mdrun.log &
echo $! > mdrun.pid
echo "mdrun started with PID $(cat mdrun.pid)"
```

**Monitor progress**:

```bash
# Check the last few log lines
tail -20 b5_qmmm.log

# Check energy so far (while still running)
echo "Total-Energy" | gmx energy -f b5_qmmm.edr -o current_energy.xvg 2>/dev/null
```

---

## 8. Checkpoint and Restart

A 500 ps run takes approximately 20 hours. If the simulation is interrupted
(network drop, power outage, job scheduler limit), you can restart from the
last checkpoint:

```bash
gmx mdrun -deffnm b5_qmmm -cpi b5_qmmm.cpt -pin off -ntmpi 1 -ntomp 8
```

**Key points**:
- Checkpoints are written at regular intervals (default: every 15 minutes;
  controllable with `gmx mdrun -cpt <minutes>`)
- By default, restart **appends** to existing output files (`.edr`, `.log`,
  `.xtc`). This means `gmx energy` sees one continuous energy file.
- Use `-noappend` if you want separate part files instead.
- NVE restart is bitwise-reproducible (no stochastic elements).

---

## 9. Analysis

### 9.1 Extract energy data

```bash
# Total energy (for NVE drift)
echo "Total-Energy" | gmx energy -f b5_qmmm.edr -o total_energy.xvg

# Temperature
echo "Temperature" | gmx energy -f b5_qmmm.edr -o temperature.xvg

# Quantum (DFTB+) energy
echo "Quantum-En." | gmx energy -f b5_qmmm.edr -o quantum_energy.xvg
```

> **Tip**: To see all available energy terms, run:
> ```bash
> echo "" | gmx energy -f b5_qmmm.edr
> ```

### 9.2 Compute NVE drift

The NVE drift measures energy conservation quality. It is computed as:

```
drift = |slope of E_tot(t)| / N_atoms    [kJ/mol/ps/atom]
```

where the slope comes from linear regression of total energy vs. time.

**Automated analysis**: Run the provided script:

```bash
bash analysis.sh b5_qmmm.edr 2652
```

This extracts all quantities and prints a pass/fail validation checklist.

### 9.3 Validation checklist

| Check | Threshold | Physical Basis |
|---|---|---|
| NVE drift | < 0.01 kJ/mol/ps/atom | specs.md S21.1; bounds accumulated energy error |
| Temperature mean | 300 +/- 10 K | NVT equilibration target |
| Temperature sigma | < 10 K | ~2x canonical fluctuation (4.7 K for N_dof=7953) |
| QM energy sigma | < 100 kJ/mol | Bounded QM potential energy surface |
| SCC anomalies | 0 frames | All SCC cycles must converge |

---

## 10. Common Pitfalls

| Problem | Symptom | Root Cause | Fix |
|---|---|---|---|
| Missing `switch-width` | NVE drift 1-10 kJ/mol/ps/atom | Hard cutoff: discontinuous forces at r_cut | Set `qmmm-cp2k-dftbplus-switch-width = 0.2` |
| Missing `damping-sigma` | SCC convergence failure; simulation crash at random time | Bare 1/r diverges at short QM-MM distance | Set `qmmm-cp2k-dftbplus-damping-sigma = 0.1` |
| `ReadInitialCharges = Yes` | Fortran ERROR STOP at startup; no error message | DFTB+ looks for `charges.bin` which doesn't exist | Set `ReadInitialCharges = No` in HSD |
| `ConvergentSCCOnly = Yes` | Random crash during MD; Fortran ERROR STOP | Single SCC failure causes unrecoverable abort | Set `ConvergentSCCOnly = No` in HSD |
| Multiple MPI ranks | Fortran ERROR STOP at startup | Each MPI rank initializes DFTB+ independently | Always use `-ntmpi 1` (QM/MM requires single rank) |
| `-pin auto` with background jobs | CPU usage drops from ~1000% to ~100% mid-run | Thread pinning collapse: all threads on 1 core | Always use `-pin off` with `gmx mdrun` |
| SK files not found | Fortran ERROR STOP at startup | DFTB+ resolves paths from CWD, not HSD location | Run `gmx mdrun` from directory with `slako/` |
| Box too small for cutoff | Unphysical energies; potential double-counting | Minimum image: L_box must be > 2 * r_cut | Ensure box edge > 2.4 nm for 1.2 nm cutoff |
| Wrong QM charge in MDP | Incorrect embedding potential | Charge mismatch between MDP and actual QM region | Set `qmmm-cp2k-qmcharge = 0` for neutral QM water |
| QM bonds constrained | Energy violation; LINCS fights DFTB+ forces | LINCS competing with quantum forces | B5: handled automatically (see Section 4, step 3) |
| Large SCC mixing parameter | SCC oscillation and divergence | Anderson/Broyden mixer overshoots | Use `MixingParameter = 0.05` or lower |

---

## 11. Expected Results

All values below come from the archived US-043b production run (500 ps NVE
with Gaussian damping). See `expected_results.md` for full provenance.

| Quantity | Value | Threshold |
|---|---|---|
| NVE drift | -4.53 x 10^-5 kJ/mol/ps/atom | < 0.01 (221x margin) |
| Total energy mean | -38,982.79 kJ/mol | — |
| Total energy sigma | 17.36 kJ/mol | — |
| Temperature | 306.26 +/- 4.77 K | 300 +/- 10 K |
| QM energy | -10,704.68 +/- 0.92 kJ/mol | sigma < 100 |
| Anomalous frames | 0 | 0 |
| Wall time | ~20 hours (500 ps) | — |

These results demonstrate excellent energy conservation: the drift is 221
times below the acceptance threshold, temperature fluctuations match the
canonical ensemble prediction, and the QM energy is remarkably stable
(0.92 kJ/mol standard deviation over 100,001 frames).

---

## 12. DFTB+ Output Files

During simulation, DFTB+ writes several files to the working directory.
With the recommended HSD settings (`WriteDetailedOut = No`, `WriteResultsTag = No`),
only these files appear:

| File | Description | Action |
|---|---|---|
| `dftb_pin.hsd` | Processed input (auto-generated) | Informational; can delete |
| `charges.bin` | Persisted Mulliken charges | Used internally; do not delete during run |
| `band.out` | Band structure info | Informational; can delete |

If you enable `WriteDetailedOut = Yes` for debugging, `detailed.out` and
`results.tag` will also appear. These are overwritten every MD step and only
reflect the last calculation.

---

## Further Reading

- Theory derivations: `docs/theory/US-044/` (parameter justification, unit analysis, numerical methods)
- Software specification: `docs/specification/specs.md` (interface contracts, module architecture)
- Benchmark validation: `docs/verification/US-043b.md` (500 ps NVE stability proof)
