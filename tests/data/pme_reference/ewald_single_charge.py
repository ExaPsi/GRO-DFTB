#!/usr/bin/env python3
"""
R-046-1: Direct Ewald summation for a single charge in a cubic box.

ThDD:T-US-046-1.4 -- Real-space Ewald sum
ThDD:T-US-046-1.5 -- Reciprocal-space Ewald sum
SDD:specs.md:§8.2  -- PME-compatible embedding

System: 1 charge Q = +1.0 e at origin in cubic box L = 3.0 nm
Probe point: R_A = (0.3, 0.4, 0.5) nm

Computes Φ_Ewald(R_A) = Φ^real + Φ^recip
(No self-correction since probe ≠ source; tin-foil boundary, no dipole correction.)

The Ewald potential from a charge Q at position R_J sensed at R_A is:
  Φ(R_A) = Q × Σ_n erfc(α|R_AJ + nL|) / |R_AJ + nL|                [real]
          + Q × (1/πV) Σ_{k≠0} (4π/k²) exp(-k²/4α²) cos(k·R_AJ)   [recip]

Units: We compute in GROMACS internal units: nm, kJ/mol.
  Potential in kJ/(mol·e) = c_one4PiEps0 × Q/r with r in nm.
  Then convert to Hartree: Φ_Ha = Φ_kJ/(mol·e) / (c_hartree2Kj × c_avogadro)

Output: ewald_single_charge.json with full provenance.
"""

import json
import numpy as np
from scipy.special import erfc

# --- GROMACS constants (CODATA 2018, from api/legacy/include/gromacs/math/units.h) ---
c_avogadro = 6.02214076e23          # /mol (exact)
c_bohr2Nm = 0.0529177210903        # nm
c_one4PiEps0 = 138.93545764498155  # kJ/(mol·nm·e²) — Coulomb constant in GROMACS units

# c_hartree2Kj in GROMACS is per molecule (not per mole):
# c_hartree2Kj = 4.3597447222071e-21 kJ
# c_hartree2Kj × c_avogadro = 2625.4996 kJ/mol
c_hartreeKjMol = 2625.4996394799    # kJ/mol per Hartree (c_hartree2Kj × c_avogadro)

# Conversion: Φ(kJ/(mol·e)) → Φ(Hartree)
pot_gmx_to_ha = 1.0 / c_hartreeKjMol

# Sanity check: unit charge at 1 Bohr should give 1 Hartree
phi_test = c_one4PiEps0 / c_bohr2Nm  # kJ/(mol·e)
phi_test_ha = phi_test * pot_gmx_to_ha
assert abs(phi_test_ha - 1.0) < 1e-6, f"Sanity check failed: phi(1 Bohr) = {phi_test_ha}, expected 1.0"

# --- System Parameters ---
L = 3.0  # Box length in nm
Q = 1.0  # Charge in elementary charges (e)
R_source = np.array([0.0, 0.0, 0.0])  # Source charge position (nm)
R_probe = np.array([0.3, 0.4, 0.5])   # Probe point (nm)

# Ewald parameters matching GROMACS defaults:
rcoulomb = 1.0   # nm
ewald_rtol = 1e-5

def calc_ewaldcoeff_q(rc, rtol):
    """Compute Ewald splitting parameter α matching GROMACS calc_ewaldcoeff_q."""
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

alpha = calc_ewaldcoeff_q(rcoulomb, ewald_rtol)  # 1/nm
print(f"Ewald splitting parameter alpha = {alpha:.10f} 1/nm")
print(f"  erfc(alpha * rcoulomb) = {erfc(alpha * rcoulomb):.6e} (should be ~ ewald_rtol)")

# --- Convergence parameters ---
nmax_real = 5     # ±5 images in each direction
kmax_recip = 20   # ±20 k-vectors in each direction

# --- Real-space sum ---
# Φ^real(R_A) = c_one4PiEps0 × Σ_n Q × erfc(α|r+nL|) / |r+nL|
# where the sum is over periodic images n, and r = R_probe - R_source
phi_real = 0.0
R_AJ = R_probe - R_source
for nx in range(-nmax_real, nmax_real + 1):
    for ny in range(-nmax_real, nmax_real + 1):
        for nz in range(-nmax_real, nmax_real + 1):
            r_vec = R_AJ + np.array([nx, ny, nz]) * L
            r = np.linalg.norm(r_vec)
            if r < 1e-15:
                continue
            phi_real += Q * erfc(alpha * r) / r

# Convert from e/nm to kJ/(mol·e)
phi_real_kjmol = c_one4PiEps0 * phi_real

print(f"Φ^real = {phi_real_kjmol:.15e} kJ/(mol·e)")

# --- Reciprocal-space sum ---
# Φ^recip(R_A) = c_one4PiEps0 × (1/V) Σ_{k≠0} (4π/k²) exp(-k²/(4α²)) × Q × cos(k·R_AJ)
V = L**3
phi_recip = 0.0
for kx in range(-kmax_recip, kmax_recip + 1):
    for ky in range(-kmax_recip, kmax_recip + 1):
        for kz in range(-kmax_recip, kmax_recip + 1):
            if kx == 0 and ky == 0 and kz == 0:
                continue
            k_vec = 2.0 * np.pi / L * np.array([kx, ky, kz])
            k2 = np.dot(k_vec, k_vec)
            structure_factor = Q * np.cos(np.dot(k_vec, R_AJ))
            phi_recip += (4.0 * np.pi / k2) * np.exp(-k2 / (4.0 * alpha**2)) * structure_factor

phi_recip /= V
# Convert from e/nm to kJ/(mol·e)
phi_recip_kjmol = c_one4PiEps0 * phi_recip

print(f"Φ^recip = {phi_recip_kjmol:.15e} kJ/(mol·e)")

# --- Total ---
phi_total_kjmol = phi_real_kjmol + phi_recip_kjmol
phi_total_ha = phi_total_kjmol * pot_gmx_to_ha
phi_real_ha = phi_real_kjmol * pot_gmx_to_ha
phi_recip_ha = phi_recip_kjmol * pot_gmx_to_ha

print(f"\nΦ_total = {phi_total_ha:.15e} Hartree")
print(f"Φ^real  = {phi_real_ha:.15e} Hartree")
print(f"Φ^recip = {phi_recip_ha:.15e} Hartree")
print(f"Φ_total = {phi_total_kjmol:.10f} kJ/(mol·e)")

# --- Electric field E = -∇Φ via finite differences ---
delta = 1e-7  # nm (central differences)
E_field_kjmol = np.zeros(3)  # kJ/(mol·e·nm)

def ewald_potential_at(R_A):
    """Compute full Ewald potential at R_A in kJ/(mol·e)."""
    R_AJ_local = R_A - R_source
    phi_r = 0.0
    for nx in range(-nmax_real, nmax_real + 1):
        for ny in range(-nmax_real, nmax_real + 1):
            for nz in range(-nmax_real, nmax_real + 1):
                r_vec = R_AJ_local + np.array([nx, ny, nz]) * L
                r = np.linalg.norm(r_vec)
                if r < 1e-15:
                    continue
                phi_r += Q * erfc(alpha * r) / r

    phi_k = 0.0
    for kx in range(-kmax_recip, kmax_recip + 1):
        for ky in range(-kmax_recip, kmax_recip + 1):
            for kz in range(-kmax_recip, kmax_recip + 1):
                if kx == 0 and ky == 0 and kz == 0:
                    continue
                k_vec = 2.0 * np.pi / L * np.array([kx, ky, kz])
                k2 = np.dot(k_vec, k_vec)
                sf = Q * np.cos(np.dot(k_vec, R_AJ_local))
                phi_k += (4.0 * np.pi / k2) * np.exp(-k2 / (4.0 * alpha**2)) * sf
    phi_k /= V
    return c_one4PiEps0 * (phi_r + phi_k)

for d in range(3):
    R_plus = R_probe.copy()
    R_minus = R_probe.copy()
    R_plus[d] += delta
    R_minus[d] -= delta
    phi_plus = ewald_potential_at(R_plus)
    phi_minus = ewald_potential_at(R_minus)
    E_field_kjmol[d] = -(phi_plus - phi_minus) / (2.0 * delta)

# Convert to Hartree/Bohr: E_Ha/Bohr = E_kJ/(mol·e·nm) × c_bohr2Nm / c_hartreeKjMol
E_field_ha_bohr = E_field_kjmol * c_bohr2Nm / c_hartreeKjMol

print(f"\nE_field (kJ/(mol·e·nm)) = [{E_field_kjmol[0]:.12e}, {E_field_kjmol[1]:.12e}, {E_field_kjmol[2]:.12e}]")
print(f"E_field (Ha/Bohr) = [{E_field_ha_bohr[0]:.12e}, {E_field_ha_bohr[1]:.12e}, {E_field_ha_bohr[2]:.12e}]")

# --- Archive results ---
results = {
    "provenance": {
        "script": "tests/data/pme_reference/ewald_single_charge.py",
        "method": "Direct Ewald summation (real + reciprocal)",
        "date": "2026-02-08",
        "python": "3.x with numpy + scipy",
        "precision": "float64",
        "gromacs_constants": "CODATA 2018 (api/legacy/include/gromacs/math/units.h)"
    },
    "system": {
        "description": "Single charge Q=+1e at origin, cubic box L=3.0 nm",
        "box_length_nm": float(L),
        "charge_e": float(Q),
        "source_position_nm": R_source.tolist(),
        "probe_position_nm": R_probe.tolist()
    },
    "parameters": {
        "rcoulomb_nm": float(rcoulomb),
        "ewald_rtol": float(ewald_rtol),
        "alpha_invnm": float(alpha),
        "nmax_real": int(nmax_real),
        "kmax_recip": int(kmax_recip)
    },
    "results": {
        "potential_hartree": float(phi_total_ha),
        "potential_real_hartree": float(phi_real_ha),
        "potential_recip_hartree": float(phi_recip_ha),
        "potential_kjmol_e": float(phi_total_kjmol),
        "potential_real_kjmol_e": float(phi_real_kjmol),
        "potential_recip_kjmol_e": float(phi_recip_kjmol),
        "electric_field_ha_bohr": E_field_ha_bohr.tolist(),
        "electric_field_kjmol_e_nm": E_field_kjmol.tolist()
    },
    "unit_conversions": {
        "c_one4PiEps0_kj_mol_nm_e2": float(c_one4PiEps0),
        "c_hartreeKjMol": float(c_hartreeKjMol),
        "c_bohr2Nm": float(c_bohr2Nm),
        "pot_gmx_to_ha": float(pot_gmx_to_ha)
    }
}

output_file = "ewald_single_charge.json"
with open(output_file, "w") as f:
    json.dump(results, f, indent=2)
print(f"\nResults archived to {output_file}")
