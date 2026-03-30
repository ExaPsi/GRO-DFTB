/*
 * SDD:specs.md:§6.1 — Unit conversion constants for GRO-DFTB
 * ThDD:CODATA2022 — All fundamental constants from NIST CODATA 2022
 *
 * Internal units: atomic units (Bohr, Hartree, e)
 * GROMACS units:  nm, kJ/mol, e
 *
 * Convert at module boundaries only; never double-convert.
 */

#ifndef GRODFTB_UNITS_H
#define GRODFTB_UNITS_H

#ifdef __cplusplus
extern "C" {
#endif

/* -----------------------------------------------------------------------
 * Length: Bohr <-> nm
 * ThDD:CODATA2022:a0 (T-US-006-1.1)
 * NIST CODATA 2022: a_0 = 0.052917721083(18) nm
 * ----------------------------------------------------------------------- */
#define GRODFTB_BOHR_TO_NM       0.052917721083
#define GRODFTB_NM_TO_BOHR       (1.0 / GRODFTB_BOHR_TO_NM)

/* -----------------------------------------------------------------------
 * Energy: Hartree <-> kJ/mol
 * ThDD:CODATA2022:Eh (T-US-006-2.1)
 * E_h = 4.3597447222060(48) aJ; E_h * N_A * 1e-3 = 2625.4996394799 kJ/mol
 * ----------------------------------------------------------------------- */
#define GRODFTB_HARTREE_TO_KJMOL  2625.4996394799
#define GRODFTB_KJMOL_TO_HARTREE  (1.0 / GRODFTB_HARTREE_TO_KJMOL)

/* -----------------------------------------------------------------------
 * Force: Hartree/Bohr <-> kJ/mol/nm
 * ThDD:derived:Eh/a0 (T-US-006-2.2)
 * ----------------------------------------------------------------------- */
#define GRODFTB_FORCE_AU_TO_GMX  (GRODFTB_HARTREE_TO_KJMOL / GRODFTB_BOHR_TO_NM)
#define GRODFTB_FORCE_GMX_TO_AU  (GRODFTB_KJMOL_TO_HARTREE * GRODFTB_BOHR_TO_NM)

/* -----------------------------------------------------------------------
 * Electrostatic potential: Hartree/e <-> kJ/mol/e
 * ----------------------------------------------------------------------- */
#define GRODFTB_POTENTIAL_AU_TO_GMX  GRODFTB_HARTREE_TO_KJMOL
#define GRODFTB_POTENTIAL_GMX_TO_AU  GRODFTB_KJMOL_TO_HARTREE

/* -----------------------------------------------------------------------
 * Electric field: Hartree/(e*Bohr) <-> kJ/(mol*e*nm)
 * ----------------------------------------------------------------------- */
#define GRODFTB_FIELD_AU_TO_GMX  (GRODFTB_HARTREE_TO_KJMOL / GRODFTB_BOHR_TO_NM)
#define GRODFTB_FIELD_GMX_TO_AU  (GRODFTB_KJMOL_TO_HARTREE * GRODFTB_BOHR_TO_NM)

/* -----------------------------------------------------------------------
 * Charge: elementary charge — no conversion needed
 * ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
 * GROMACS dual-constant convention note (SDD:specs.md:§6.1)
 *
 * GROMACS uses CODATA 2018 constants (c_bohr2Nm = 0.0529177210903).
 * GRO-DFTB uses CODATA 2022 internally. At the GROMACS interface
 * boundary (Module 3), GROMACS constants must be used to avoid
 * double-rounding. This header's constants are for libgrodftb internals.
 * ----------------------------------------------------------------------- */
#define GRODFTB_GROMACS_CODATA_NOTE \
    "GROMACS uses CODATA 2018 constants (c_bohr2Nm = 0.0529177210903). " \
    "GRO-DFTB uses CODATA 2022 internally. At the GROMACS interface, "   \
    "use GROMACS constants to avoid double-rounding."

/* -----------------------------------------------------------------------
 * In-place conversion functions (SDD:specs.md:§6.2)
 * ----------------------------------------------------------------------- */

/** Convert coordinate array in-place: nm -> Bohr. */
void grodftb_coords_nm_to_bohr(int n, double *coords);

/** Convert coordinate array in-place: Bohr -> nm. */
void grodftb_coords_bohr_to_nm(int n, double *coords);

/** Convert force array in-place: Hartree/Bohr -> kJ/mol/nm. */
void grodftb_forces_au_to_gmx(int n, double *forces);

/** Convert force array in-place: kJ/mol/nm -> Hartree/Bohr. */
void grodftb_forces_gmx_to_au(int n, double *forces);

#ifdef __cplusplus
}
#endif

#endif /* GRODFTB_UNITS_H */
