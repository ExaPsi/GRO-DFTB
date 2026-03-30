/*
 * SDD:specs.md:§9 — Link Atom Handler (Module 5)
 * ThDD:T-US-033-2.1 — Link atom placement theory
 *
 * Handle covalent bonds crossing the QM/MM boundary by placing hydrogen
 * cap atoms ("link atoms") and projecting forces to real atoms.
 *
 * LGPL-3.0-or-later
 */

#ifndef GRODFTB_LINKATOM_H
#define GRODFTB_LINKATOM_H

#include "grodftb/error.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---------------------------------------------------------------------------
 * SDD:specs.md:§9.4 — Charge redistribution schemes
 * --------------------------------------------------------------------------- */
#define GRODFTB_CHARGE_NONE   0  /* No charge modification */
#define GRODFTB_CHARGE_ZERO   1  /* Set MM boundary atom charge to zero */
#define GRODFTB_CHARGE_SHIFT  2  /* Redistribute charge to MM neighbors */

/* ---------------------------------------------------------------------------
 * US-035: Additional error codes for charge redistribution
 * --------------------------------------------------------------------------- */
#define GRODFTB_ERR_SINGULAR_MATRIX   20  /* G matrix is singular */
#define GRODFTB_ERR_NO_M2_NEIGHBORS   21  /* No M2 neighbors for shift scheme */

/* ---------------------------------------------------------------------------
 * US-035: Numerical thresholds for charge redistribution
 * ThDD:T-US-035-N.2 — Singularity detection
 * --------------------------------------------------------------------------- */
#define GRODFTB_SINGULAR_THRESHOLD    1e-12  /* |det(G)| threshold */
#define GRODFTB_CONDITION_THRESHOLD   1e12   /* Condition number limit */

/* ---------------------------------------------------------------------------
 * SDD:specs.md:§9.5 — Link atom definition
 *
 * Stores the mapping from a QM boundary atom to its MM partner, plus
 * the parameters needed for link atom placement and force projection.
 * --------------------------------------------------------------------------- */
typedef struct {
    int    qm_atom;       /* Index of QM atom A in QM-only array (0-based) */
    int    mm_atom;       /* Global index of MM atom B */
    int    link_species;  /* Species index for link atom (H), 0-based */
    double g_factor;      /* Current scaling factor g = d_link / d_AB */
    double ref_bond_len;  /* Reference A-H bond length (Bohr) */
} grodftb_link_atom_t;

/* ---------------------------------------------------------------------------
 * Opaque handle for the link atom handler
 * --------------------------------------------------------------------------- */
typedef struct grodftb_linkatom_handler *grodftb_linkatom_handle_t;

/* ---------------------------------------------------------------------------
 * Creation and destruction
 * --------------------------------------------------------------------------- */

/**
 * Create a link atom handler from frontier definitions.
 *
 * SDD:specs.md:§9.6 — Link atoms identified from topology (done by caller)
 *
 * @param nlinks        Number of link atoms (may be 0)
 * @param qm_atoms      QM atom indices for each link (0-based in QM array) [nlinks]
 * @param mm_atoms      MM atom global indices for each link [nlinks]
 * @param ref_lengths   Reference A-H bond lengths in Bohr [nlinks]
 * @param h_species_idx Species index for hydrogen in DFTB+ HSD
 * @param charge_scheme GRODFTB_CHARGE_NONE/ZERO/SHIFT
 * @param handle_out    Receives the new handle on success
 * @return GRODFTB_SUCCESS on success, error code otherwise
 */
int grodftb_linkatom_create(int nlinks,
                            const int *qm_atoms,
                            const int *mm_atoms,
                            const double *ref_lengths,
                            int h_species_idx,
                            int charge_scheme,
                            grodftb_linkatom_handle_t *handle_out);

/**
 * Release all resources associated with the link atom handler.
 *
 * @param handle  Pointer to handle (set to NULL on return)
 */
void grodftb_linkatom_destroy(grodftb_linkatom_handle_t *handle);

/**
 * Get number of link atoms.
 *
 * @param handle  Link atom handler (may be NULL)
 * @return Number of link atoms, or 0 if handle is NULL
 */
int grodftb_linkatom_count(grodftb_linkatom_handle_t handle);

/**
 * Get the current g-factor for a specific link atom.
 *
 * ThDD:T-US-033-4.4 — g = d_link / d_AB
 *
 * Note: g-factor is updated by compute_positions(). Call that first.
 *
 * @param handle    Link atom handler
 * @param link_idx  Link atom index (0-based)
 * @param g_out     Receives the g-factor
 * @return GRODFTB_SUCCESS on success
 */
int grodftb_linkatom_get_g_factor(grodftb_linkatom_handle_t handle,
                                  int link_idx,
                                  double *g_out);

/* ---------------------------------------------------------------------------
 * Link atom position computation
 * --------------------------------------------------------------------------- */

/**
 * Compute link atom positions from current QM and MM coordinates.
 *
 * ThDD:T-US-033-4.1 — R_L = R_A + d_link * unit_vec(R_B - R_A)
 * ThDD:T-US-033-4.4 — g = d_link / d_AB (updated internally)
 *
 * All coordinates are in atomic units (Bohr).
 *
 * @param handle        Link atom handler
 * @param nqm           Number of real QM atoms
 * @param qm_coords     QM coordinates in Bohr, flat [3*nqm]
 * @param mm_coords     MM coordinates in Bohr, flat array (global indexing)
 * @param link_pos_out  Receives link atom positions in Bohr, flat [3*nlinks]
 * @return GRODFTB_SUCCESS on success, GRODFTB_ERR_INVALID_ARGUMENT if d_AB < epsilon
 */
int grodftb_linkatom_compute_positions(grodftb_linkatom_handle_t handle,
                                       int nqm,
                                       const double *qm_coords,
                                       const double *mm_coords,
                                       double *link_pos_out);

/**
 * Build augmented coordinate array: real QM atoms + link atoms.
 *
 * Output array layout: [qm0, qm1, ..., qm_{nqm-1}, link0, link1, ...]
 * where each entry is 3 consecutive doubles (x, y, z).
 *
 * @param handle           Link atom handler
 * @param nqm              Number of real QM atoms
 * @param qm_coords        QM coordinates in Bohr, flat [3*nqm]
 * @param mm_coords        MM coordinates in Bohr, flat array (global indexing)
 * @param augmented_out    Pre-allocated output [3*(nqm + nlinks)]
 * @return GRODFTB_SUCCESS on success
 */
int grodftb_linkatom_augment_coords(grodftb_linkatom_handle_t handle,
                                    int nqm,
                                    const double *qm_coords,
                                    const double *mm_coords,
                                    double *augmented_out);

/**
 * Build augmented species array: real QM species + link atom species.
 *
 * Output layout: [sp0, sp1, ..., sp_{nqm-1}, h_sp, h_sp, ...]
 *
 * @param handle           Link atom handler
 * @param nqm              Number of real QM atoms
 * @param qm_species       QM species indices [nqm]
 * @param augmented_out    Pre-allocated output [nqm + nlinks]
 * @return GRODFTB_SUCCESS on success
 */
int grodftb_linkatom_augment_species(grodftb_linkatom_handle_t handle,
                                     int nqm,
                                     const int *qm_species,
                                     int *augmented_out);

/* ---------------------------------------------------------------------------
 * Force projection
 * --------------------------------------------------------------------------- */

/**
 * Project link atom forces back to real atoms A and B.
 *
 * ThDD:T-US-033-5.9 — GROMACS force spreading formula:
 *   F_B = g * (F_L - (F_L . r_hat) * r_hat)
 *   F_A = F_L - F_B
 *
 * This correctly handles the chain-rule terms from variable g.
 *
 * Input: forces on all atoms including links, layout [F_QM..., F_L...]
 * The QM forces from positions 0..3*nqm-1 are copied to output.
 * Link atom forces are projected and added to the appropriate QM and MM atoms.
 *
 * All forces are in atomic units (Hartree/Bohr).
 *
 * @param handle           Link atom handler (g_factors must be current)
 * @param nqm              Number of real QM atoms
 * @param qm_coords        QM coordinates in Bohr [3*nqm] (for projection calculation)
 * @param mm_coords        MM coordinates in Bohr (global indexing)
 * @param augmented_forces Forces in Hartree/Bohr [3*(nqm + nlinks)]
 * @param qm_forces_out    QM forces output, Hartree/Bohr [3*nqm] (overwritten, then link contributions added)
 * @param mm_forces_out    MM forces output, Hartree/Bohr (global indexing, link contributions added)
 * @return GRODFTB_SUCCESS on success
 */
int grodftb_linkatom_project_forces(grodftb_linkatom_handle_t handle,
                                    int nqm,
                                    const double *qm_coords,
                                    const double *mm_coords,
                                    const double *augmented_forces,
                                    double *qm_forces_out,
                                    double *mm_forces_out);

/* ---------------------------------------------------------------------------
 * US-035: Charge Redistribution
 * SDD:specs.md:§9.4 — Charge redistribution schemes for link atom boundaries
 * ThDD:T-US-035-2.1 through T-US-035-6.13
 * --------------------------------------------------------------------------- */

/**
 * US-035: Result structure for charge redistribution.
 * SDD:specs.md:§9.4
 *
 * Contains the modified charges, force corrections (for shift scheme),
 * and diagnostic information about the redistribution.
 */
typedef struct grodftb_redistrib_result {
    double *modified_charges;   /**< Output charges [ncharges] (caller must not free) */
    double *force_corrections;  /**< Chain-rule forces (shift only) [3*ncharges] (caller must not free) */
    int     fallback_used;      /**< 1 if reduced constraints used (N<4 or singular) */
    double  charge_error;       /**< |sum(Q') - sum(Q)| for verification */
    double  dipole_error[3];    /**< |sum(dQ*d) - 0| per component (shift scheme) */
    /* Internal fields */
    int     ncharges;           /**< Number of charges (for memory management) */
} grodftb_redistrib_result_t;

/**
 * Create a charge redistribution result structure.
 *
 * @param result_out  Receives pointer to new result structure
 * @param ncharges    Number of MM charges to handle
 * @return GRODFTB_SUCCESS on success, GRODFTB_ERR_ALLOC_FAILED on memory error
 */
int grodftb_redistrib_result_create(grodftb_redistrib_result_t **result_out,
                                    int ncharges);

/**
 * Destroy a charge redistribution result structure.
 *
 * @param result  Pointer to result structure (set to NULL on return)
 */
void grodftb_redistrib_result_destroy(grodftb_redistrib_result_t **result);

/**
 * Set M2 neighbor information for a link atom.
 *
 * This function specifies which MM atoms are neighbors of the M1 atom
 * (boundary MM atom) for each link. Required for the shift scheme.
 *
 * @param handle       Link atom handler
 * @param link_idx     Link atom index (0-based)
 * @param n_neighbors  Number of M2 neighbors for this link
 * @param neighbors    Array of global MM indices for M2 atoms [n_neighbors]
 * @return GRODFTB_SUCCESS on success
 */
int grodftb_linkatom_set_m2_neighbors(grodftb_linkatom_handle_t handle,
                                       int link_idx,
                                       int n_neighbors,
                                       const int *neighbors);

/**
 * Get the number of M2 neighbors for a link atom.
 *
 * @param handle     Link atom handler
 * @param link_idx   Link atom index (0-based)
 * @param count_out  Receives the neighbor count
 * @return GRODFTB_SUCCESS on success
 */
int grodftb_linkatom_get_m2_neighbor_count(grodftb_linkatom_handle_t handle,
                                           int link_idx,
                                           int *count_out);

/**
 * Redistribute MM charges near QM/MM boundary.
 *
 * ThDD:T-US-035-2.1 (none), T-US-035-3.1 (zero), T-US-035-4.1-4.22 (shift)
 *
 * Applies the charge redistribution scheme configured in the handler.
 * The result structure contains:
 * - modified_charges: The charges after redistribution
 * - force_corrections: Chain-rule force corrections (shift scheme only)
 * - fallback_used: Flag indicating if reduced constraints were used
 * - charge_error: Total charge conservation error (should be < 1e-12)
 * - dipole_error: Dipole moment error per component (shift scheme)
 *
 * All quantities are in atomic units (Bohr for positions, e for charges,
 * Hartree/Bohr for forces).
 *
 * @param handle       Link atom handler with M2 neighbors configured
 * @param ncharges     Number of MM charges
 * @param charges      Original MM charges [ncharges]
 * @param charge_pos   MM charge positions in Bohr [3*ncharges]
 * @param qm_charges   Mulliken charges on QM atoms [n_qm]
 * @param qm_pos       QM atom positions in Bohr [3*n_qm]
 * @param n_qm         Number of QM atoms
 * @param result       Pre-allocated result structure (from grodftb_redistrib_result_create)
 * @return GRODFTB_SUCCESS on success, error code otherwise
 */
int grodftb_linkatom_redistribute_charges(grodftb_linkatom_handle_t handle,
                                          int ncharges,
                                          const double *charges,
                                          const double *charge_pos,
                                          const double *qm_charges,
                                          const double *qm_pos,
                                          int n_qm,
                                          grodftb_redistrib_result_t *result);

/**
 * Get the charge redistribution scheme from a handler.
 *
 * @param handle      Link atom handler
 * @param scheme_out  Receives GRODFTB_CHARGE_NONE/ZERO/SHIFT
 * @return GRODFTB_SUCCESS on success
 */
int grodftb_linkatom_get_charge_scheme(grodftb_linkatom_handle_t handle,
                                       int *scheme_out);

#ifdef __cplusplus
}
#endif

#endif /* GRODFTB_LINKATOM_H */
