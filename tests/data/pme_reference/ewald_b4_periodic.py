#!/usr/bin/env python3
"""
R-046-2: Direct Ewald summation for B4 water dimer + 1 MM charge in periodic box.

ThDD:T-US-046-1.4 -- Real-space Ewald sum
ThDD:T-US-046-1.5 -- Reciprocal-space Ewald sum
SDD:specs.md:§8.2  -- PME-compatible embedding

System: B4 water dimer (6 QM atoms) + 1 MM charge (Q = +1.0 e) in cubic box L = 2.0 nm.
QM atom coordinates from tests/data/b4/geo.gen (Bohr → nm).
MM charge at (-3.0, 0.0, 0.0) Angstrom = (-0.3, 0.0, 0.0) nm.

Computes Φ_Ewald(R_A) at all 6 QM atom sites from the single MM charge.

Output: ewald_b4_periodic.json with full provenance.
"""

import json
import numpy as np
from scipy.special import erfc

# --- GROMACS constants (CODATA 2018) ---
c_avogadro = 6.02214076e23
c_bohr2Nm = 0.0529177210903
c_one4PiEps0 = 138.93545764498155  # kJ/(mol·nm·e²)
c_hartreeKjMol = 2625.4996394799    # kJ/mol per Hartree

pot_gmx_to_ha = 1.0 / c_hartreeKjMol

# --- System Parameters ---
L = 2.0  # Box length in nm

# QM atoms from tests/data/b4/geo.gen (Bohr coordinates)
qm_atoms_bohr = np.array([
    [-0.702196, -0.056060, 0.009942],   # O (water 1)
    [-1.022193, 0.846776, -0.011489],    # H (water 1)
    [0.257521, 0.042121, 0.005219],      # H (water 1)
    [2.220871, 0.026717, 0.000620],      # O (water 2)
    [2.597493, -0.411663, 0.766745],     # H (water 2)
    [2.593135, -0.449496, -0.744782],    # H (water 2)
])
qm_labels = ["O1", "H1a", "H1b", "O2", "H2a", "H2b"]
qm_atoms_nm = qm_atoms_bohr * c_bohr2Nm
nQM = len(qm_atoms_nm)

# MM charge at (-3.0, 0.0, 0.0) Angstrom = (-0.3, 0.0, 0.0) nm
mm_charges = [{"position_nm": [-0.3, 0.0, 0.0], "charge_e": 1.0}]

# Ewald parameters
rcoulomb = 0.9   # nm (typical for smaller box)
ewald_rtol = 1e-5

def calc_ewaldcoeff_q(rc, rtol):
    """Compute Ewald splitting parameter α matching GROMACS."""
    low = 0.0
    high = 10.0 / rc
    for _ in range(100):
        beta = 0.5 * (low + high)
        if erfc(beta * rc) > rtol:
            low = beta
        else:
            high = beta
        if high - low < 1e-14:
            break
    return beta

alpha = calc_ewaldcoeff_q(rcoulomb, ewald_rtol)
print(f"Ewald splitting parameter alpha = {alpha:.10f} 1/nm")

# Convergence parameters
nmax_real = 5
kmax_recip = 25  # More k-vectors for smaller box

V = L**3

def ewald_potential_and_field(R_A, mm_charges_list, L, alpha, nmax_real, kmax_recip):
    """Compute full Ewald potential and field at R_A from MM charges."""
    V = L**3
    phi_real = 0.0
    phi_recip = 0.0
    E_real = np.zeros(3)
    E_recip = np.zeros(3)

    for mm in mm_charges_list:
        R_J = np.array(mm["position_nm"])
        Q_J = mm["charge_e"]
        R_AJ = R_A - R_J

        # Real-space sum
        for nx in range(-nmax_real, nmax_real + 1):
            for ny in range(-nmax_real, nmax_real + 1):
                for nz in range(-nmax_real, nmax_real + 1):
                    r_vec = R_AJ + np.array([nx, ny, nz]) * L
                    r = np.linalg.norm(r_vec)
                    if r < 1e-15:
                        continue
                    erfc_val = erfc(alpha * r)
                    phi_real += Q_J * erfc_val / r
                    # Field: E = -∇Φ = Q_J × K(r,α) × r_hat
                    # K(r,α) = erfc(αr)/r³ + 2α/√π exp(-α²r²)/r²
                    K = erfc_val / r**3 + 2.0 * alpha / np.sqrt(np.pi) * np.exp(-alpha**2 * r**2) / r**2
                    E_real += Q_J * K * r_vec  # This is actually -dΦ/dR = +Q*K*r_hat, pointing away

        # Reciprocal-space sum
        for kx in range(-kmax_recip, kmax_recip + 1):
            for ky in range(-kmax_recip, kmax_recip + 1):
                for kz in range(-kmax_recip, kmax_recip + 1):
                    if kx == 0 and ky == 0 and kz == 0:
                        continue
                    k_vec = 2.0 * np.pi / L * np.array([kx, ky, kz])
                    k2 = np.dot(k_vec, k_vec)
                    prefac = (4.0 * np.pi / k2) * np.exp(-k2 / (4.0 * alpha**2))
                    kr = np.dot(k_vec, R_AJ)
                    phi_recip += Q_J * prefac * np.cos(kr)
                    # Field: dΦ/dR_A = -(Q_J/V) × (4π/k²) exp(-k²/4α²) × k × sin(k·R_AJ)
                    E_recip += Q_J * prefac * k_vec * np.sin(kr)

        phi_recip /= V
        E_recip /= V

    # Convert to GROMACS units: kJ/(mol·e)
    phi_real_kjmol = c_one4PiEps0 * phi_real
    phi_recip_kjmol = c_one4PiEps0 * phi_recip

    # Field in kJ/(mol·e·nm): E_gmx = c_one4PiEps0 × E_reduced
    # Note: E_real is the gradient of 1/r terms (not negative gradient),
    # so the electric field E = -∇Φ has the OPPOSITE sign convention.
    # Actually: E_real = +Q_J × K(r,α) × r_hat points AWAY from source for +Q.
    # This IS the correct E = -∇Φ direction (potential decreases away from +Q).
    # Wait, for Φ = Q/r, ∇Φ = -Q*r_hat/r², so E = -∇Φ = +Q*r_hat/r² (points away from +Q).
    # For erfc: ∇(erfc(αr)/r) = -(erfc(αr)/r³ + 2α/√π exp(-α²r²)/r²) × r_vec
    # So E = -∇Φ = +Q_J × K(r,α) × r_vec (divided by r for hat). But K already includes /r³ etc.
    # Let me redo: Φ = Q_J × erfc(αr)/r
    # dΦ/dx = Q_J × d(erfc(αr)/r)/dx = Q_J × [-erfc(αr)/r³ - 2α/√π exp(-α²r²)/r²] × x
    # So E_x = -dΦ/dx = Q_J × [erfc(αr)/r³ + 2α/√π exp(-α²r²)/r²] × x = Q_J × K × x
    # This is positive when x > 0 and Q > 0 (field points away from positive charge). Correct.

    E_real_kjmol = c_one4PiEps0 * E_real
    E_recip_kjmol = c_one4PiEps0 * E_recip

    return phi_real_kjmol, phi_recip_kjmol, E_real_kjmol, E_recip_kjmol

# --- Compute at all QM sites ---
print(f"\n{'Atom':<6} {'Φ_real (Ha)':<22} {'Φ_recip (Ha)':<22} {'Φ_total (Ha)':<22}")
print("-" * 74)

qm_results = []
for i in range(nQM):
    phi_r, phi_k, E_r, E_k = ewald_potential_and_field(
        qm_atoms_nm[i], mm_charges, L, alpha, nmax_real, kmax_recip
    )

    phi_total = phi_r + phi_k
    E_total = E_r + E_k

    phi_r_ha = phi_r * pot_gmx_to_ha
    phi_k_ha = phi_k * pot_gmx_to_ha
    phi_total_ha = phi_total * pot_gmx_to_ha

    E_total_ha_bohr = E_total * c_bohr2Nm / c_hartreeKjMol

    print(f"{qm_labels[i]:<6} {phi_r_ha:<22.15e} {phi_k_ha:<22.15e} {phi_total_ha:<22.15e}")

    qm_results.append({
        "atom_index": i,
        "atom_label": qm_labels[i],
        "position_nm": qm_atoms_nm[i].tolist(),
        "position_bohr": qm_atoms_bohr[i].tolist(),
        "potential_total_hartree": float(phi_total_ha),
        "potential_real_hartree": float(phi_r_ha),
        "potential_recip_hartree": float(phi_k_ha),
        "potential_total_kjmol_e": float(phi_total),
        "electric_field_ha_bohr": E_total_ha_bohr.tolist(),
        "electric_field_kjmol_e_nm": E_total.tolist()
    })

# --- Also compute FD-validated field at atom 0 ---
print(f"\nFD field validation at atom 0 ({qm_labels[0]}):")
delta = 1e-7  # nm
E_fd = np.zeros(3)
for d in range(3):
    R_plus = qm_atoms_nm[0].copy()
    R_minus = qm_atoms_nm[0].copy()
    R_plus[d] += delta
    R_minus[d] -= delta
    phi_r_p, phi_k_p, _, _ = ewald_potential_and_field(R_plus, mm_charges, L, alpha, nmax_real, kmax_recip)
    phi_r_m, phi_k_m, _, _ = ewald_potential_and_field(R_minus, mm_charges, L, alpha, nmax_real, kmax_recip)
    E_fd[d] = -((phi_r_p + phi_k_p) - (phi_r_m + phi_k_m)) / (2.0 * delta)

phi_r_0, phi_k_0, E_r_0, E_k_0 = ewald_potential_and_field(
    qm_atoms_nm[0], mm_charges, L, alpha, nmax_real, kmax_recip
)
E_analytic = E_r_0 + E_k_0
print(f"  Analytic: E = [{E_analytic[0]:.12e}, {E_analytic[1]:.12e}, {E_analytic[2]:.12e}] kJ/(mol·e·nm)")
print(f"  FD:       E = [{E_fd[0]:.12e}, {E_fd[1]:.12e}, {E_fd[2]:.12e}] kJ/(mol·e·nm)")
rel_err = np.abs(E_analytic - E_fd) / (np.abs(E_analytic) + 1e-30)
print(f"  Rel err:  [{rel_err[0]:.6e}, {rel_err[1]:.6e}, {rel_err[2]:.6e}]")

# --- Archive ---
results = {
    "provenance": {
        "script": "tests/data/pme_reference/ewald_b4_periodic.py",
        "method": "Direct Ewald summation (real + reciprocal)",
        "date": "2026-02-08",
        "python": "3.x with numpy + scipy",
        "precision": "float64",
        "gromacs_constants": "CODATA 2018",
        "qm_geometry_source": "tests/data/b4/geo.gen (Bohr)",
        "mm_charge_source": "tests/data/b4/dftb_in.hsd PointCharges"
    },
    "system": {
        "description": "B4 water dimer + 1 MM charge in cubic box L=2.0 nm",
        "box_length_nm": float(L),
        "nQM": int(nQM),
        "qm_atom_labels": qm_labels,
        "qm_atoms_bohr": qm_atoms_bohr.tolist(),
        "qm_atoms_nm": qm_atoms_nm.tolist(),
        "mm_charges": mm_charges
    },
    "parameters": {
        "rcoulomb_nm": float(rcoulomb),
        "ewald_rtol": float(ewald_rtol),
        "alpha_invnm": float(alpha),
        "nmax_real": int(nmax_real),
        "kmax_recip": int(kmax_recip)
    },
    "results": qm_results
}

output_file = "ewald_b4_periodic.json"
with open(output_file, "w") as f:
    json.dump(results, f, indent=2)
print(f"\nResults archived to {output_file}")
