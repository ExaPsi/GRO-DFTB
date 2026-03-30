#!/usr/bin/env python3
"""
Compute Delta_PBC: the difference between full Ewald and minimum-image
bare Coulomb potentials at QM sites from MM charges.

For the double-counting correction in the charge-inserted Ewald approach,
we use minimum-image bare Coulomb 1/r. The actual periodic potential is
the Ewald sum. The difference is Delta_PBC.

Usage: python3 compute_delta_pbc.py <gro_file>
"""
import numpy as np
import json
import sys
from scipy.special import erfc

# Constants
NM_TO_BOHR = 1.0 / 0.0529177210903
HARTREE_TO_KJMOL = 2625.4996394799

def read_gro(filename):
    """Read GRO file, return coordinates in nm and box vectors."""
    with open(filename) as f:
        lines = f.readlines()
    natom = int(lines[1].strip())
    coords = np.zeros((natom, 3))
    for i in range(natom):
        line = lines[2 + i]
        coords[i, 0] = float(line[20:28])
        coords[i, 1] = float(line[28:36])
        coords[i, 2] = float(line[36:44])
    box_line = lines[2 + natom].split()
    box = np.array([float(x) for x in box_line[:3]])
    return coords, box

def minimum_image_distance(r1, r2, box):
    """Compute minimum-image distance between two points."""
    dr = r2 - r1
    dr -= box * np.round(dr / box)
    return np.sqrt(np.sum(dr**2)), dr

def compute_bare_coulomb(qm_coords, mm_coords, mm_charges, box):
    """Compute bare Coulomb potential at QM sites from MM charges using minimum-image.
    Returns potential in Hartree (atomic units: q/r with r in Bohr)."""
    n_qm = len(qm_coords)
    n_mm = len(mm_coords)
    pot = np.zeros(n_qm)

    for i in range(n_qm):
        for j in range(n_mm):
            dist_nm, _ = minimum_image_distance(qm_coords[i], mm_coords[j], box)
            dist_bohr = dist_nm * NM_TO_BOHR
            if dist_bohr > 1e-6:
                pot[i] += mm_charges[j] / dist_bohr  # Hartree

    return pot

def compute_ewald_potential(qm_coords, mm_coords, mm_charges, box, alpha=None, kmax=10):
    """
    Compute Ewald sum potential at QM sites from MM charges.

    alpha: Ewald splitting parameter (in 1/nm). If None, use GROMACS default.
    kmax: max reciprocal vector index.
    """
    if alpha is None:
        # For rcoulomb=1.0 nm and ewald_rtol=1e-5:  alpha ≈ 3.12 nm^-1
        alpha = 3.123413  # nm^-1

    n_qm = len(qm_coords)
    n_mm = len(mm_coords)
    V = box[0] * box[1] * box[2]  # nm^3

    # Real-space sum (with minimum-image convention, which is exact for Ewald
    # when alpha is chosen so that erfc(alpha * L/2) << 1)
    # Accumulate in e/nm, then convert to Hartree
    pot_real = np.zeros(n_qm)
    for i in range(n_qm):
        for j in range(n_mm):
            dist, _ = minimum_image_distance(qm_coords[i], mm_coords[j], box)
            if dist > 1e-10:
                pot_real[i] += mm_charges[j] * erfc(alpha * dist) / dist  # e/nm

    # Convert from e/nm to Hartree:
    # Φ(Ha) = q(e) / r(Bohr) = q(e) / (r(nm) * NM_TO_BOHR)
    # So (q/r in e/nm) / NM_TO_BOHR = Φ in Hartree
    pot_real /= NM_TO_BOHR

    # Reciprocal-space sum (in e/nm, then convert)
    pot_recip = np.zeros(n_qm)
    twopi = 2.0 * np.pi

    for nx in range(-kmax, kmax + 1):
        for ny in range(-kmax, kmax + 1):
            for nz in range(-kmax, kmax + 1):
                if nx == 0 and ny == 0 and nz == 0:
                    continue
                k = twopi * np.array([nx / box[0], ny / box[1], nz / box[2]])  # nm^-1
                k2 = np.sum(k**2)

                # Structure factor from MM charges
                # S(k) = sum_j q_j * exp(-i k.r_j)
                S_real = 0.0
                S_imag = 0.0
                for j in range(n_mm):
                    phase = np.dot(k, mm_coords[j])
                    S_real += mm_charges[j] * np.cos(phase)
                    S_imag += mm_charges[j] * np.sin(phase)

                # Potential at QM sites (in e/nm)
                prefactor = (4.0 * np.pi / V) * np.exp(-k2 / (4.0 * alpha**2)) / k2
                for i in range(n_qm):
                    phase_i = np.dot(k, qm_coords[i])
                    pot_recip[i] += prefactor * (S_real * np.cos(phase_i) + S_imag * np.sin(phase_i))

    # Convert from e/nm to Hartree
    pot_recip /= NM_TO_BOHR

    # No self-term needed (QM and MM are different sets)
    # No dipole correction for tin-foil boundary conditions

    return pot_real, pot_recip

def main():
    gro_file = sys.argv[1] if len(sys.argv) > 1 else \
        "tests/data/b5_pme_drift/initial.gro"

    print(f"Reading: {gro_file}")
    coords, box = read_gro(gro_file)
    print(f"Box: {box[0]:.4f} x {box[1]:.4f} x {box[2]:.4f} nm")
    print(f"Total atoms: {len(coords)}")

    # QM atoms: indices 1182, 1183, 1184 (0-based) = atoms 1183-1185 (1-based)
    qm_indices = np.array([1182, 1183, 1184])

    # MM atoms: all others
    all_indices = np.arange(len(coords))
    mm_indices = np.setdiff1d(all_indices, qm_indices)

    qm_coords = coords[qm_indices]
    mm_coords = coords[mm_indices]

    # TIP3P charges (QM charges zeroed by preprocessing)
    # MM charges: OW=-0.834, HW1=+0.417, HW2=+0.417 for each water
    mm_charges = np.zeros(len(mm_indices))
    for i, idx in enumerate(mm_indices):
        atom_in_mol = idx % 3  # 0=OW, 1=HW1, 2=HW2
        if atom_in_mol == 0:
            mm_charges[i] = -0.834
        else:
            mm_charges[i] = 0.417

    print(f"\nQM atoms: {qm_indices + 1} (1-based)")
    print(f"MM atoms: {len(mm_indices)} ({len(mm_indices)//3} waters)")
    print(f"Total MM charge: {mm_charges.sum():.6f} e")

    # Compute bare Coulomb (minimum-image)
    print("\nComputing bare Coulomb (minimum-image)...")
    pot_bare = compute_bare_coulomb(qm_coords, mm_coords, mm_charges, box)

    # Compute Ewald sum
    print("Computing Ewald sum (alpha=3.123413/nm, kmax=10)...")
    pot_ewald_real, pot_ewald_recip = compute_ewald_potential(
        qm_coords, mm_coords, mm_charges, box, alpha=3.123413, kmax=10
    )
    pot_ewald = pot_ewald_real + pot_ewald_recip

    # Delta_PBC = Ewald - bare
    delta_pbc = pot_ewald - pot_bare

    labels = ["O (QM)", "H1 (QM)", "H2 (QM)"]

    print("\n" + "=" * 70)
    print(f"{'Atom':>10s} {'Φ_bare (Ha)':>14s} {'Φ_Ewald (Ha)':>14s} {'Δ_PBC (Ha)':>14s} {'Δ_PBC (kJ/mol)':>16s}")
    print("=" * 70)
    for i in range(3):
        print(f"{labels[i]:>10s} {pot_bare[i]:14.6f} {pot_ewald[i]:14.6f} "
              f"{delta_pbc[i]:14.6f} {delta_pbc[i]*HARTREE_TO_KJMOL:16.4f}")

    print(f"\n{'Mean |Δ_PBC|':>10s} {np.mean(np.abs(delta_pbc)):14.6f} Ha = "
          f"{np.mean(np.abs(delta_pbc)) * HARTREE_TO_KJMOL:.4f} kJ/mol")

    print(f"\nEwald decomposition:")
    print(f"  Real-space:      {pot_ewald_real.mean():14.6f} Ha (mean)")
    print(f"  Reciprocal:      {pot_ewald_recip.mean():14.6f} Ha (mean)")
    print(f"  Total Ewald:     {pot_ewald.mean():14.6f} Ha (mean)")
    print(f"  Bare Coulomb:    {pot_bare.mean():14.6f} Ha (mean)")

    # Also compute the correction energy from Delta_PBC
    # Using typical QM Mulliken charges for water: O ~ -0.6e, H ~ +0.3e
    qm_charges_typical = np.array([-0.6, 0.3, 0.3])
    E_delta = np.sum(qm_charges_typical * delta_pbc)
    print(f"\n  Correction energy from Δ_PBC (using q_O=-0.6, q_H=+0.3):")
    print(f"    ΔE_PBC = Σ q_A × Δ_PBC_A = {E_delta:.6f} Ha = {E_delta*HARTREE_TO_KJMOL:.4f} kJ/mol")

    # Relative to bare potential
    print(f"\n  Relative magnitude: |Δ_PBC| / |Φ_bare| = {np.mean(np.abs(delta_pbc))/np.mean(np.abs(pot_bare)):.4f}")

    # Save results
    results = {
        "gro_file": gro_file,
        "box_nm": box.tolist(),
        "alpha_inv_nm": 3.123413,
        "kmax": 10,
        "qm_atoms_1based": (qm_indices + 1).tolist(),
        "n_mm_atoms": int(len(mm_indices)),
        "phi_bare_Ha": pot_bare.tolist(),
        "phi_ewald_real_Ha": pot_ewald_real.tolist(),
        "phi_ewald_recip_Ha": pot_ewald_recip.tolist(),
        "phi_ewald_total_Ha": pot_ewald.tolist(),
        "delta_pbc_Ha": delta_pbc.tolist(),
        "delta_pbc_mean_abs_Ha": float(np.mean(np.abs(delta_pbc))),
        "delta_pbc_mean_abs_kjmol": float(np.mean(np.abs(delta_pbc)) * HARTREE_TO_KJMOL),
        "relative_magnitude": float(np.mean(np.abs(delta_pbc))/np.mean(np.abs(pot_bare))),
        "energy_correction_Ha": float(E_delta),
        "energy_correction_kjmol": float(E_delta * HARTREE_TO_KJMOL),
    }

    outfile = gro_file.replace('.gro', '_delta_pbc.json')
    if 'initial' in outfile:
        outfile = "docs/session/20260213/delta_pbc_results.json"
    with open(outfile, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {outfile}")


if __name__ == "__main__":
    main()
