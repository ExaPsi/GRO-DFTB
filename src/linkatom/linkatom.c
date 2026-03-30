/*
 * SDD:specs.md:§9 — Link Atom Handler (Module 5) Implementation
 * ThDD:T-US-033-2.1, T-US-033-4.1, T-US-033-5.9 — Link atom placement and force projection
 * ThDD:T-US-035-2.1 through T-US-035-6.13 — Charge redistribution schemes
 *
 * Handles covalent bonds crossing the QM/MM boundary by:
 * 1. Placing hydrogen cap atoms along the QM-MM bond vector
 * 2. Projecting link atom forces to the real boundary atoms
 * 3. Redistributing M1 charges to M2 neighbors (shift scheme)
 *
 * Units: All internal computations in atomic units (Bohr, Hartree).
 *        Convert at module boundaries only.
 *
 * LGPL-3.0-or-later
 */

#include "grodftb/linkatom.h"
#include "grodftb/error.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Minimum distance threshold to avoid division by zero.
 * ThDD:T-US-033-N1 — d_AB > epsilon, else unphysical configuration.
 * Value: ~10^-10 nm in Bohr = ~2e-9 Bohr (machine epsilon scale) */
#define LINKATOM_MIN_DISTANCE 1e-10

/* ---------------------------------------------------------------------------
 * Internal data structure for M2 neighbor tracking (US-035)
 * --------------------------------------------------------------------------- */
typedef struct {
    int   n_neighbors;    /* Number of M2 neighbors for this link */
    int  *neighbors;      /* Global MM indices of M2 atoms [n_neighbors] */
} grodftb_m2_info_t;

/* ---------------------------------------------------------------------------
 * Internal data structure
 * --------------------------------------------------------------------------- */
struct grodftb_linkatom_handler {
    int                  nlinks;        /* Number of link atoms */
    grodftb_link_atom_t *links;         /* Array of link atom definitions */
    int                  charge_scheme; /* GRODFTB_CHARGE_NONE/ZERO/SHIFT */
    int                  h_species;     /* Species index for hydrogen */
    grodftb_m2_info_t   *m2_info;       /* M2 neighbor info per link [nlinks], US-035 */
};

/* ---------------------------------------------------------------------------
 * Helper: 3D vector operations (inline for performance)
 * --------------------------------------------------------------------------- */

/* Compute v = b - a */
static inline void vec3_sub(const double *a, const double *b, double *v)
{
    v[0] = b[0] - a[0];
    v[1] = b[1] - a[1];
    v[2] = b[2] - a[2];
}

/* Compute Euclidean norm */
static inline double vec3_norm(const double *v)
{
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

/* Compute dot product */
static inline double vec3_dot(const double *a, const double *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

/* Scale vector: v *= s */
static inline void vec3_scale(double *v, double s)
{
    v[0] *= s;
    v[1] *= s;
    v[2] *= s;
}

/* Add scaled vector: a += s * b */
static inline void vec3_add_scaled(double *a, const double *b, double s)
{
    a[0] += s * b[0];
    a[1] += s * b[1];
    a[2] += s * b[2];
}

/* ---------------------------------------------------------------------------
 * Creation and destruction
 * --------------------------------------------------------------------------- */

int grodftb_linkatom_create(int nlinks,
                            const int *qm_atoms,
                            const int *mm_atoms,
                            const double *ref_lengths,
                            int h_species_idx,
                            int charge_scheme,
                            grodftb_linkatom_handle_t *handle_out)
{
    /* SDD:specs.md:§18.1 — Validate arguments */
    if (handle_out == NULL) {
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    *handle_out = NULL;

    /* Zero links is valid (no covalent boundary) */
    if (nlinks < 0) {
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }

    if (nlinks > 0 && (qm_atoms == NULL || mm_atoms == NULL || ref_lengths == NULL)) {
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }

    /* Allocate handler */
    struct grodftb_linkatom_handler *h = malloc(sizeof(*h));
    if (h == NULL) {
        return GRODFTB_ERR_ALLOC_FAILED;
    }

    h->nlinks = nlinks;
    h->charge_scheme = charge_scheme;
    h->h_species = h_species_idx;
    h->links = NULL;
    h->m2_info = NULL;

    if (nlinks > 0) {
        h->links = malloc(nlinks * sizeof(grodftb_link_atom_t));
        if (h->links == NULL) {
            free(h);
            return GRODFTB_ERR_ALLOC_FAILED;
        }

        /* Allocate M2 neighbor info array (US-035) */
        h->m2_info = malloc(nlinks * sizeof(grodftb_m2_info_t));
        if (h->m2_info == NULL) {
            free(h->links);
            free(h);
            return GRODFTB_ERR_ALLOC_FAILED;
        }

        /* Initialize link atom definitions and M2 info */
        for (int i = 0; i < nlinks; i++) {
            h->links[i].qm_atom = qm_atoms[i];
            h->links[i].mm_atom = mm_atoms[i];
            h->links[i].link_species = h_species_idx;
            h->links[i].ref_bond_len = ref_lengths[i];
            h->links[i].g_factor = 0.0;  /* Will be set by compute_positions */

            /* Initialize M2 info to empty (US-035) */
            h->m2_info[i].n_neighbors = 0;
            h->m2_info[i].neighbors = NULL;
        }
    }

    *handle_out = h;
    return GRODFTB_SUCCESS;
}

void grodftb_linkatom_destroy(grodftb_linkatom_handle_t *handle)
{
    if (handle == NULL || *handle == NULL) {
        return;
    }

    struct grodftb_linkatom_handler *h = *handle;

    /* Free M2 neighbor arrays (US-035) */
    if (h->m2_info != NULL) {
        for (int i = 0; i < h->nlinks; i++) {
            if (h->m2_info[i].neighbors != NULL) {
                free(h->m2_info[i].neighbors);
            }
        }
        free(h->m2_info);
    }

    if (h->links != NULL) {
        free(h->links);
    }
    free(h);
    *handle = NULL;
}

int grodftb_linkatom_count(grodftb_linkatom_handle_t handle)
{
    if (handle == NULL) {
        return 0;
    }
    return handle->nlinks;
}

int grodftb_linkatom_get_g_factor(grodftb_linkatom_handle_t handle,
                                  int link_idx,
                                  double *g_out)
{
    if (handle == NULL || g_out == NULL) {
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    if (link_idx < 0 || link_idx >= handle->nlinks) {
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }

    *g_out = handle->links[link_idx].g_factor;
    return GRODFTB_SUCCESS;
}

/* ---------------------------------------------------------------------------
 * Link atom position computation
 *
 * ThDD:T-US-033-4.1 — R_L = R_A + d_link * unit_vec(R_B - R_A)
 *
 * Algorithm LA-1 from 03_numerical_methods.md:
 * 1. Compute displacement vector: r = R_B - R_A
 * 2. Compute distance: d_AB = ||r||
 * 3. Compute unit vector: u = r / d_AB
 * 4. Compute link position: R_L = R_A + d_link * u
 * --------------------------------------------------------------------------- */

int grodftb_linkatom_compute_positions(grodftb_linkatom_handle_t handle,
                                       int nqm,
                                       const double *qm_coords,
                                       const double *mm_coords,
                                       double *link_pos_out)
{
    /* nqm reserved for future bounds checking; silence unused warning */
    (void)nqm;

    if (handle == NULL) {
        return GRODFTB_ERR_NULL_HANDLE;
    }
    if (handle->nlinks == 0) {
        return GRODFTB_SUCCESS;  /* No work to do */
    }
    if (qm_coords == NULL || mm_coords == NULL || link_pos_out == NULL) {
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }

    for (int i = 0; i < handle->nlinks; i++) {
        grodftb_link_atom_t *link = &handle->links[i];

        /* Get positions of boundary atoms */
        const double *r_A = &qm_coords[3 * link->qm_atom];
        const double *r_B = &mm_coords[3 * link->mm_atom];

        /* Step 1: Displacement vector r = R_B - R_A */
        double r[3];
        vec3_sub(r_A, r_B, r);

        /* Step 2: Distance d_AB = ||r|| */
        double d_AB = vec3_norm(r);

        /*
         * ThDD:T-US-033-N1 — Check for collapsed atoms.
         * This is unphysical and would cause division by zero.
         */
        if (d_AB < LINKATOM_MIN_DISTANCE) {
            return GRODFTB_ERR_INVALID_ARGUMENT;
        }

        /* Step 3: Unit vector and g-factor
         * ThDD:T-US-033-4.4 — g = d_link / d_AB */
        double inv_d_AB = 1.0 / d_AB;
        double g = link->ref_bond_len * inv_d_AB;

        /* Store g-factor for force projection */
        link->g_factor = g;

        /* Step 4: Link position R_L = R_A + d_link * u = R_A + g * r
         * ThDD:T-US-033-4.1 */
        double *r_L = &link_pos_out[3 * i];
        r_L[0] = r_A[0] + g * r[0];
        r_L[1] = r_A[1] + g * r[1];
        r_L[2] = r_A[2] + g * r[2];
    }

    return GRODFTB_SUCCESS;
}

int grodftb_linkatom_augment_coords(grodftb_linkatom_handle_t handle,
                                    int nqm,
                                    const double *qm_coords,
                                    const double *mm_coords,
                                    double *augmented_out)
{
    if (handle == NULL) {
        return GRODFTB_ERR_NULL_HANDLE;
    }
    if (qm_coords == NULL || augmented_out == NULL) {
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    if (handle->nlinks > 0 && mm_coords == NULL) {
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }

    /* Copy QM coordinates to output */
    memcpy(augmented_out, qm_coords, 3 * nqm * sizeof(double));

    /* Compute link atom positions directly into the output array */
    if (handle->nlinks > 0) {
        double *link_start = &augmented_out[3 * nqm];
        int rc = grodftb_linkatom_compute_positions(handle, nqm, qm_coords,
                                                    mm_coords, link_start);
        if (rc != GRODFTB_SUCCESS) {
            return rc;
        }
    }

    return GRODFTB_SUCCESS;
}

int grodftb_linkatom_augment_species(grodftb_linkatom_handle_t handle,
                                     int nqm,
                                     const int *qm_species,
                                     int *augmented_out)
{
    if (handle == NULL) {
        return GRODFTB_ERR_NULL_HANDLE;
    }
    if (qm_species == NULL || augmented_out == NULL) {
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }

    /* Copy QM species to output */
    memcpy(augmented_out, qm_species, nqm * sizeof(int));

    /* Append H species for all link atoms */
    for (int i = 0; i < handle->nlinks; i++) {
        augmented_out[nqm + i] = handle->h_species;
    }

    return GRODFTB_SUCCESS;
}

/* ---------------------------------------------------------------------------
 * Force projection
 *
 * ThDD:T-US-033-5.9 — GROMACS force spreading formula:
 *   F_B = g * (F_L - (F_L . r_hat) * r_hat)
 *   F_A = F_L - F_B
 *
 * This is the energy-consistent projection that correctly handles the
 * chain-rule terms from the variable g-factor.
 *
 * Algorithm LA-2 from 03_numerical_methods.md:
 * 1. Compute displacement with PBC: r = pbc_dx(R_B, R_A)  [no PBC here yet]
 * 2. Compute distance: d_AB = ||r||
 * 3. Compute inverse distance: inv_d = 1 / d_AB
 * 4. Compute g factor: g = d_link * inv_d
 * 5. Compute projected force component: F_parallel = (F_L . r) * inv_d^2
 * 6. Compute force on MM atom: F_B = g * (F_L - F_parallel * r)
 * 7. Compute force on QM atom: F_A = F_L - F_B
 * --------------------------------------------------------------------------- */

int grodftb_linkatom_project_forces(grodftb_linkatom_handle_t handle,
                                    int nqm,
                                    const double *qm_coords,
                                    const double *mm_coords,
                                    const double *augmented_forces,
                                    double *qm_forces_out,
                                    double *mm_forces_out)
{
    if (handle == NULL) {
        return GRODFTB_ERR_NULL_HANDLE;
    }
    if (augmented_forces == NULL || qm_forces_out == NULL) {
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    if (handle->nlinks > 0 && (qm_coords == NULL || mm_coords == NULL || mm_forces_out == NULL)) {
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }

    /* Copy QM forces from augmented array */
    memcpy(qm_forces_out, augmented_forces, 3 * nqm * sizeof(double));

    /* Project link atom forces */
    for (int i = 0; i < handle->nlinks; i++) {
        const grodftb_link_atom_t *link = &handle->links[i];

        /* Get link atom force */
        const double *F_L = &augmented_forces[3 * (nqm + i)];

        /* Get boundary atom positions */
        const double *r_A = &qm_coords[3 * link->qm_atom];
        const double *r_B = &mm_coords[3 * link->mm_atom];

        /* Step 1: Displacement r = R_B - R_A */
        double r[3];
        vec3_sub(r_A, r_B, r);

        /* Step 2-3: Distance and inverse */
        double d_AB = vec3_norm(r);
        if (d_AB < LINKATOM_MIN_DISTANCE) {
            /* Should not happen if compute_positions was called first */
            return GRODFTB_ERR_INVALID_ARGUMENT;
        }
        double inv_d = 1.0 / d_AB;

        /* Step 4: g factor (use stored value for consistency) */
        double g = link->g_factor;

        /* Step 5: Projected force component
         * F_parallel = (F_L . r) * inv_d^2 */
        double F_dot_r = vec3_dot(F_L, r);
        double proj_coeff = F_dot_r * inv_d * inv_d;

        /* Step 6: Force on MM atom B
         * ThDD:T-US-033-5.9 — F_B = g * (F_L - (F_L . r_hat) * r_hat)
         *                        = g * (F_L - proj_coeff * r) */
        double F_B[3];
        F_B[0] = g * (F_L[0] - proj_coeff * r[0]);
        F_B[1] = g * (F_L[1] - proj_coeff * r[1]);
        F_B[2] = g * (F_L[2] - proj_coeff * r[2]);

        /* Step 7: Force on QM atom A
         * ThDD:T-US-033-5.9 — F_A = F_L - F_B */
        double F_A[3];
        F_A[0] = F_L[0] - F_B[0];
        F_A[1] = F_L[1] - F_B[1];
        F_A[2] = F_L[2] - F_B[2];

        /* Accumulate to output arrays */
        double *qm_F_A = &qm_forces_out[3 * link->qm_atom];
        qm_F_A[0] += F_A[0];
        qm_F_A[1] += F_A[1];
        qm_F_A[2] += F_A[2];

        double *mm_F_B = &mm_forces_out[3 * link->mm_atom];
        mm_F_B[0] += F_B[0];
        mm_F_B[1] += F_B[1];
        mm_F_B[2] += F_B[2];
    }

    return GRODFTB_SUCCESS;
}

/* ===========================================================================
 * US-035: Charge Redistribution Implementation
 * ThDD:T-US-035-2.1 through T-US-035-6.13
 * ===========================================================================*/

/* ---------------------------------------------------------------------------
 * Accessor for charge scheme
 * --------------------------------------------------------------------------- */

int grodftb_linkatom_get_charge_scheme(grodftb_linkatom_handle_t handle,
                                       int *scheme_out)
{
    if (handle == NULL || scheme_out == NULL) {
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    *scheme_out = handle->charge_scheme;
    return GRODFTB_SUCCESS;
}

/* ---------------------------------------------------------------------------
 * M2 neighbor management
 * --------------------------------------------------------------------------- */

int grodftb_linkatom_set_m2_neighbors(grodftb_linkatom_handle_t handle,
                                       int link_idx,
                                       int n_neighbors,
                                       const int *neighbors)
{
    if (handle == NULL) {
        return GRODFTB_ERR_NULL_HANDLE;
    }
    if (link_idx < 0 || link_idx >= handle->nlinks) {
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    if (n_neighbors < 0) {
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    if (n_neighbors > 0 && neighbors == NULL) {
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }

    grodftb_m2_info_t *info = &handle->m2_info[link_idx];

    /* Free existing neighbors if any */
    if (info->neighbors != NULL) {
        free(info->neighbors);
        info->neighbors = NULL;
    }
    info->n_neighbors = 0;

    /* Allocate and copy new neighbors */
    if (n_neighbors > 0) {
        info->neighbors = malloc(n_neighbors * sizeof(int));
        if (info->neighbors == NULL) {
            return GRODFTB_ERR_ALLOC_FAILED;
        }
        memcpy(info->neighbors, neighbors, n_neighbors * sizeof(int));
        info->n_neighbors = n_neighbors;
    }

    return GRODFTB_SUCCESS;
}

int grodftb_linkatom_get_m2_neighbor_count(grodftb_linkatom_handle_t handle,
                                           int link_idx,
                                           int *count_out)
{
    if (handle == NULL || count_out == NULL) {
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    if (link_idx < 0 || link_idx >= handle->nlinks) {
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }

    *count_out = handle->m2_info[link_idx].n_neighbors;
    return GRODFTB_SUCCESS;
}

/* ---------------------------------------------------------------------------
 * Result structure management
 * --------------------------------------------------------------------------- */

int grodftb_redistrib_result_create(grodftb_redistrib_result_t **result_out,
                                    int ncharges)
{
    if (result_out == NULL) {
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    *result_out = NULL;

    if (ncharges < 0) {
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }

    grodftb_redistrib_result_t *r = malloc(sizeof(*r));
    if (r == NULL) {
        return GRODFTB_ERR_ALLOC_FAILED;
    }

    r->ncharges = ncharges;
    r->modified_charges = NULL;
    r->force_corrections = NULL;
    r->fallback_used = 0;
    r->charge_error = 0.0;
    r->dipole_error[0] = 0.0;
    r->dipole_error[1] = 0.0;
    r->dipole_error[2] = 0.0;

    if (ncharges > 0) {
        r->modified_charges = malloc(ncharges * sizeof(double));
        if (r->modified_charges == NULL) {
            free(r);
            return GRODFTB_ERR_ALLOC_FAILED;
        }

        r->force_corrections = malloc(3 * ncharges * sizeof(double));
        if (r->force_corrections == NULL) {
            free(r->modified_charges);
            free(r);
            return GRODFTB_ERR_ALLOC_FAILED;
        }

        /* Zero-initialize */
        memset(r->modified_charges, 0, ncharges * sizeof(double));
        memset(r->force_corrections, 0, 3 * ncharges * sizeof(double));
    }

    *result_out = r;
    return GRODFTB_SUCCESS;
}

void grodftb_redistrib_result_destroy(grodftb_redistrib_result_t **result)
{
    if (result == NULL || *result == NULL) {
        return;
    }

    grodftb_redistrib_result_t *r = *result;
    if (r->modified_charges != NULL) {
        free(r->modified_charges);
    }
    if (r->force_corrections != NULL) {
        free(r->force_corrections);
    }
    free(r);
    *result = NULL;
}

/* ---------------------------------------------------------------------------
 * Internal: Apply none scheme (passthrough)
 * ThDD:T-US-035-2.1 — Q'_M1 = Q_M1
 * --------------------------------------------------------------------------- */
static void apply_none_scheme(int ncharges,
                              const double *charges,
                              grodftb_redistrib_result_t *result)
{
    /* ThDD:T-US-035-2.1 — Direct passthrough of all charges */
    memcpy(result->modified_charges, charges, ncharges * sizeof(double));

    /* No force corrections for none scheme */
    memset(result->force_corrections, 0, 3 * ncharges * sizeof(double));

    result->fallback_used = 0;
    result->charge_error = 0.0;
    result->dipole_error[0] = 0.0;
    result->dipole_error[1] = 0.0;
    result->dipole_error[2] = 0.0;
}

/* ---------------------------------------------------------------------------
 * Internal: Apply zero scheme
 * ThDD:T-US-035-3.1 — Q'_M1 = 0
 * --------------------------------------------------------------------------- */
static void apply_zero_scheme(int ncharges,
                              const double *charges,
                              int nlinks,
                              const grodftb_link_atom_t *links,
                              grodftb_redistrib_result_t *result)
{
    /* Start with copy of original charges */
    memcpy(result->modified_charges, charges, ncharges * sizeof(double));

    /* ThDD:T-US-035-3.1 — Zero the M1 charges */
    double total_removed = 0.0;
    for (int i = 0; i < nlinks; i++) {
        int m1_idx = links[i].mm_atom;
        if (m1_idx >= 0 && m1_idx < ncharges) {
            total_removed += result->modified_charges[m1_idx];
            result->modified_charges[m1_idx] = 0.0;
        }
    }

    /* No force corrections for zero scheme (ThDD:T-US-035-3.6) */
    memset(result->force_corrections, 0, 3 * ncharges * sizeof(double));

    result->fallback_used = 0;
    /* ThDD:T-US-035-3.7 — Total charge changes by -Q_M1 */
    result->charge_error = fabs(total_removed);
    result->dipole_error[0] = 0.0;
    result->dipole_error[1] = 0.0;
    result->dipole_error[2] = 0.0;
}

/* ---------------------------------------------------------------------------
 * Internal: Build the G = A*A^T matrix (4x4)
 * ThDD:T-US-035-4.20
 * --------------------------------------------------------------------------- */
static void build_G_matrix(int n_neighbors,
                           const double *displacements,  /* [3*n_neighbors] */
                           double G[4][4])
{
    /*
     * ThDD:T-US-035-4.20 — G matrix construction:
     * G = | N      S_x    S_y    S_z   |
     *     | S_x    S_xx   S_xy   S_xz  |
     *     | S_y    S_xy   S_yy   S_yz  |
     *     | S_z    S_xz   S_yz   S_zz  |
     *
     * where S_x = sum(d_{i,x}), S_xx = sum(d_{i,x}^2), etc.
     */

    /* Initialize sums */
    double S[3] = {0.0, 0.0, 0.0};           /* S_x, S_y, S_z */
    double SS[3][3] = {{0.0}};               /* S_xx, S_xy, ... */

    /* Accumulate sums */
    for (int i = 0; i < n_neighbors; i++) {
        const double *d = &displacements[3 * i];
        for (int a = 0; a < 3; a++) {
            S[a] += d[a];
            for (int b = 0; b < 3; b++) {
                SS[a][b] += d[a] * d[b];
            }
        }
    }

    /* Build G matrix */
    G[0][0] = (double)n_neighbors;
    for (int a = 0; a < 3; a++) {
        G[0][a + 1] = S[a];
        G[a + 1][0] = S[a];
        for (int b = 0; b < 3; b++) {
            G[a + 1][b + 1] = SS[a][b];
        }
    }
}

/* ---------------------------------------------------------------------------
 * Internal: 4x4 determinant using cofactor expansion
 * ThDD:T-US-035-N.2 — For singularity detection
 * --------------------------------------------------------------------------- */
static double determinant_3x3(double m[3][3])
{
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
         - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
         + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

static double determinant_4x4(double G[4][4])
{
    double det = 0.0;

    /* Expand along first row */
    for (int j = 0; j < 4; j++) {
        /* Build 3x3 minor */
        double minor[3][3];
        for (int row = 1; row < 4; row++) {
            int col_out = 0;
            for (int col = 0; col < 4; col++) {
                if (col != j) {
                    minor[row - 1][col_out] = G[row][col];
                    col_out++;
                }
            }
        }

        double sign = (j % 2 == 0) ? 1.0 : -1.0;
        det += sign * G[0][j] * determinant_3x3(minor);
    }

    return det;
}

/* ---------------------------------------------------------------------------
 * Internal: Solve 4x4 linear system G*lambda = b using Cramer's rule
 * ThDD:T-US-035-4.21 — lambda = G^{-1} b
 *
 * Returns 0 on success, nonzero if singular.
 * --------------------------------------------------------------------------- */
static int solve_4x4_system(double G[4][4], double b[4], double lambda[4])
{
    double det_G = determinant_4x4(G);

    /* ThDD:T-US-035-N.2 — Singularity check */
    if (fabs(det_G) < GRODFTB_SINGULAR_THRESHOLD) {
        return 1;  /* Singular */
    }

    double inv_det = 1.0 / det_G;

    /* Cramer's rule: lambda_j = det(G with column j replaced by b) / det(G) */
    for (int j = 0; j < 4; j++) {
        /* Build modified matrix */
        double G_mod[4][4];
        for (int row = 0; row < 4; row++) {
            for (int col = 0; col < 4; col++) {
                G_mod[row][col] = (col == j) ? b[row] : G[row][col];
            }
        }
        lambda[j] = determinant_4x4(G_mod) * inv_det;
    }

    return 0;
}

/* ---------------------------------------------------------------------------
 * Internal: Solve 2x2 reduced constraint system (for N<4 or singular G)
 * ThDD:T-US-035-4.15a,b — Preserve charge and bond-direction dipole only
 *
 * Constraints:
 *   sum(Delta Q) = Q_M1
 *   sum(Delta Q * (d_i . n_hat)) = 0
 *
 * where n_hat = (R_A - R_M1) / |R_A - R_M1| is the bond direction.
 *
 * Minimum-norm solution with 2 constraints.
 *
 * Returns the bond-direction Lagrange multiplier (lambda2[1]) via lambda_bond_out
 * for computing chain-rule force corrections. The force formula is:
 *   d(delta_q[i])/d(R_M2_j)_mu = delta_{ij} * lambda_bond * n_hat_mu
 *   d(delta_q[i])/d(R_M1)_mu = -lambda_bond * n_hat_mu
 * --------------------------------------------------------------------------- */
static int solve_reduced_constraints(int n_neighbors,
                                     const double *displacements,
                                     double q_m1,
                                     const double *bond_dir,  /* unit vector, 3 elements */
                                     double *delta_q,         /* output [n_neighbors] */
                                     double *lambda_bond_out) /* output: bond-dir Lagrange multiplier */
{
    /* Initialize output */
    *lambda_bond_out = 0.0;

    if (n_neighbors == 0) {
        return 1;  /* Cannot satisfy any constraints */
    }

    if (n_neighbors == 1) {
        /* ThDD:T-US-035-4.18 — Single neighbor: all charge goes there */
        delta_q[0] = q_m1;
        /* No bond-direction constraint with 1 neighbor, so no lambda_bond */
        *lambda_bond_out = 0.0;
        return 0;
    }

    /*
     * Build 2xN constraint matrix A:
     * A = | 1      1      ...  1     |
     *     | d1.n   d2.n   ...  dN.n  |
     *
     * where di.n is the projection of displacement i onto bond direction.
     */
    double *row0 = malloc(n_neighbors * sizeof(double));  /* All 1s */
    double *row1 = malloc(n_neighbors * sizeof(double));  /* Projections */
    if (row0 == NULL || row1 == NULL) {
        if (row0) free(row0);
        if (row1) free(row1);
        return 2;  /* Allocation failure */
    }

    for (int i = 0; i < n_neighbors; i++) {
        row0[i] = 1.0;
        row1[i] = vec3_dot(&displacements[3 * i], bond_dir);
    }

    /*
     * G = A * A^T is 2x2:
     * G = | sum(1)       sum(di.n)     |
     *     | sum(di.n)    sum((di.n)^2) |
     */
    double G2[2][2];
    G2[0][0] = 0.0;
    G2[0][1] = 0.0;
    G2[1][0] = 0.0;
    G2[1][1] = 0.0;

    for (int i = 0; i < n_neighbors; i++) {
        G2[0][0] += row0[i] * row0[i];
        G2[0][1] += row0[i] * row1[i];
        G2[1][0] += row1[i] * row0[i];
        G2[1][1] += row1[i] * row1[i];
    }

    /* b = (Q_M1, 0)^T */
    double b2[2] = {q_m1, 0.0};

    /* Solve G2 * lambda2 = b2 */
    double det = G2[0][0] * G2[1][1] - G2[0][1] * G2[1][0];
    if (fabs(det) < GRODFTB_SINGULAR_THRESHOLD) {
        /* Even reduced system is singular - use equal distribution */
        double equal_share = q_m1 / n_neighbors;
        for (int i = 0; i < n_neighbors; i++) {
            delta_q[i] = equal_share;
        }
        /* No bond-direction constraint when singular */
        *lambda_bond_out = 0.0;
        free(row0);
        free(row1);
        return 0;
    }

    double inv_det = 1.0 / det;
    double lambda2[2];
    lambda2[0] = inv_det * (G2[1][1] * b2[0] - G2[0][1] * b2[1]);
    lambda2[1] = inv_det * (-G2[1][0] * b2[0] + G2[0][0] * b2[1]);

    /* Return the bond-direction Lagrange multiplier for force corrections */
    *lambda_bond_out = lambda2[1];

    /* Delta Q_i = A^T * lambda = lambda[0] * 1 + lambda[1] * (di.n) */
    for (int i = 0; i < n_neighbors; i++) {
        delta_q[i] = lambda2[0] * row0[i] + lambda2[1] * row1[i];
    }

    free(row0);
    free(row1);
    return 0;
}

/* ---------------------------------------------------------------------------
 * Internal: Compute QM electrostatic potential at a point
 * ThDD:T-US-035-6.5 — Phi_QM(R) = sum_alpha q_alpha / |R_alpha - R|
 * --------------------------------------------------------------------------- */
static double compute_qm_potential(int n_qm,
                                   const double *qm_charges,
                                   const double *qm_pos,
                                   const double *point)
{
    double phi = 0.0;

    for (int alpha = 0; alpha < n_qm; alpha++) {
        double r[3];
        r[0] = qm_pos[3 * alpha + 0] - point[0];
        r[1] = qm_pos[3 * alpha + 1] - point[1];
        r[2] = qm_pos[3 * alpha + 2] - point[2];

        double dist = vec3_norm(r);
        if (dist > LINKATOM_MIN_DISTANCE) {
            phi += qm_charges[alpha] / dist;
        }
    }

    return phi;
}

/* ---------------------------------------------------------------------------
 * Internal: Apply shift scheme for a single link atom
 * ThDD:T-US-035-4.1-4.22 — Full charge redistribution with dipole preservation
 * ThDD:T-US-035-6.9, 6.10 — Chain-rule force corrections
 *
 * Parameters:
 *   link        - Link atom definition
 *   m2_info     - M2 neighbor information
 *   ncharges    - Total number of MM charges
 *   charge_pos  - All MM charge positions [3*ncharges]
 *   qm_pos      - QM atom positions [3*n_qm]
 *   n_qm        - Number of QM atoms
 *   qm_charges  - QM Mulliken charges [n_qm]
 *   result      - Result structure (charges and forces updated)
 *   fallback_used - Set to 1 if reduced constraints were used
 * --------------------------------------------------------------------------- */
static int apply_shift_single_link(const grodftb_link_atom_t *link,
                                   const grodftb_m2_info_t *m2_info,
                                   int ncharges,
                                   const double *charge_pos,
                                   const double *qm_pos,
                                   int n_qm,
                                   const double *qm_charges,
                                   grodftb_redistrib_result_t *result,
                                   int *fallback_used)
{
    /* ncharges reserved for future bounds checking */
    (void)ncharges;

    int m1_idx = link->mm_atom;
    int qm_idx = link->qm_atom;
    int n_neighbors = m2_info->n_neighbors;

    /* Get M1 position and charge */
    const double *r_m1 = &charge_pos[3 * m1_idx];
    double q_m1 = result->modified_charges[m1_idx];

    /* ThDD:T-US-035-4.1 — Zero the M1 charge */
    result->modified_charges[m1_idx] = 0.0;

    /* Edge case: N = 0 (no M2 neighbors) */
    if (n_neighbors == 0) {
        /* ThDD:T-US-035-4 — Fall back to zero scheme for this link */
        *fallback_used = 1;
        return GRODFTB_SUCCESS;
    }

    /* Compute displacement vectors: d_i = R_M2_i - R_M1 */
    /* ThDD:T-US-035-4.6 */
    double *displacements = malloc(3 * n_neighbors * sizeof(double));
    double *delta_q = malloc(n_neighbors * sizeof(double));
    if (displacements == NULL || delta_q == NULL) {
        if (displacements) free(displacements);
        if (delta_q) free(delta_q);
        return GRODFTB_ERR_ALLOC_FAILED;
    }

    for (int i = 0; i < n_neighbors; i++) {
        int m2_idx = m2_info->neighbors[i];
        const double *r_m2 = &charge_pos[3 * m2_idx];
        displacements[3 * i + 0] = r_m2[0] - r_m1[0];
        displacements[3 * i + 1] = r_m2[1] - r_m1[1];
        displacements[3 * i + 2] = r_m2[2] - r_m1[2];
    }

    /* Lambda vector for force corrections (full 4-constraint case) */
    double lambda[4] = {0.0, 0.0, 0.0, 0.0};
    /* For reduced constraints, we use bond_dir and lambda_bond instead */
    double bond_dir[3] = {0.0, 0.0, 0.0};
    double lambda_bond = 0.0;
    int use_reduced_forces = 0;

    if (n_neighbors >= 4) {
        /* Full 4-constraint system */
        /* ThDD:T-US-035-4.20 — Build G matrix */
        double G[4][4];
        build_G_matrix(n_neighbors, displacements, G);

        /* ThDD:T-US-035-4.11 — Constraint vector b */
        double b[4] = {q_m1, 0.0, 0.0, 0.0};

        /* ThDD:T-US-035-4.21 — Solve for lambda */
        int singular = solve_4x4_system(G, b, lambda);

        if (singular) {
            /* Fall back to reduced constraints */
            *fallback_used = 1;
            use_reduced_forces = 1;

            /* Compute bond direction for reduced constraints */
            const double *r_qm = &qm_pos[3 * qm_idx];
            bond_dir[0] = r_qm[0] - r_m1[0];
            bond_dir[1] = r_qm[1] - r_m1[1];
            bond_dir[2] = r_qm[2] - r_m1[2];
            double bond_len = vec3_norm(bond_dir);
            if (bond_len > LINKATOM_MIN_DISTANCE) {
                bond_dir[0] /= bond_len;
                bond_dir[1] /= bond_len;
                bond_dir[2] /= bond_len;
            }

            int rc = solve_reduced_constraints(n_neighbors, displacements,
                                               q_m1, bond_dir, delta_q, &lambda_bond);
            if (rc != 0) {
                free(displacements);
                free(delta_q);
                return GRODFTB_ERR_SINGULAR_MATRIX;
            }
        } else {
            /* ThDD:T-US-035-4.22 — Compute charge shifts */
            for (int i = 0; i < n_neighbors; i++) {
                delta_q[i] = lambda[0]
                           + lambda[1] * displacements[3 * i + 0]
                           + lambda[2] * displacements[3 * i + 1]
                           + lambda[3] * displacements[3 * i + 2];
            }
        }
    } else {
        /* N < 4: Use reduced constraints */
        /* ThDD:T-US-035-4.15a,b */
        *fallback_used = 1;
        use_reduced_forces = 1;

        /* Compute bond direction */
        const double *r_qm = &qm_pos[3 * qm_idx];
        bond_dir[0] = r_qm[0] - r_m1[0];
        bond_dir[1] = r_qm[1] - r_m1[1];
        bond_dir[2] = r_qm[2] - r_m1[2];
        double bond_len = vec3_norm(bond_dir);
        if (bond_len > LINKATOM_MIN_DISTANCE) {
            bond_dir[0] /= bond_len;
            bond_dir[1] /= bond_len;
            bond_dir[2] /= bond_len;
        }

        int rc = solve_reduced_constraints(n_neighbors, displacements,
                                           q_m1, bond_dir, delta_q, &lambda_bond);
        if (rc != 0) {
            free(displacements);
            free(delta_q);
            return (rc == 2) ? GRODFTB_ERR_ALLOC_FAILED : GRODFTB_ERR_SINGULAR_MATRIX;
        }
    }

    /* Apply charge shifts to M2 atoms */
    /* ThDD:T-US-035-4.2 — Q'_M2_i = Q_M2_i + Delta Q_M2_i */
    for (int i = 0; i < n_neighbors; i++) {
        int m2_idx = m2_info->neighbors[i];
        result->modified_charges[m2_idx] += delta_q[i];
    }

    /* Compute chain-rule force corrections */
    /* ThDD:T-US-035-6.9, 6.10 */

    if (use_reduced_forces) {
        /*
         * For reduced constraints (N < 4 or singular G), the force formula is different.
         * The charge shift is: delta_q[i] = lambda_0 + lambda_bond * (d_i . n_hat)
         * where lambda_0 is the charge constraint multiplier and lambda_bond is the
         * bond-direction dipole multiplier.
         *
         * The derivatives are:
         *   d(delta_q[i])/d(R_M2_j)_mu = delta_{ij} * lambda_bond * n_hat_mu
         *   d(delta_q[i])/d(R_M1)_mu = -lambda_bond * n_hat_mu
         *
         * So the effective "lambda_spatial" for force corrections is lambda_bond * n_hat.
         */
        double lambda_spatial[3];
        lambda_spatial[0] = lambda_bond * bond_dir[0];
        lambda_spatial[1] = lambda_bond * bond_dir[1];
        lambda_spatial[2] = lambda_bond * bond_dir[2];

        /* Sum of QM potential at M2 sites for M1 force */
        double sum_phi_lambda[3] = {0.0, 0.0, 0.0};

        for (int i = 0; i < n_neighbors; i++) {
            int m2_idx = m2_info->neighbors[i];
            const double *r_m2 = &charge_pos[3 * m2_idx];

            /* ThDD:T-US-035-6.5 — QM potential at M2 site */
            double phi_qm = compute_qm_potential(n_qm, qm_charges, qm_pos, r_m2);

            /* ThDD:T-US-035-6.9 — Force on M2: F_M2^(chain) = -Phi_QM * lambda_spatial */
            result->force_corrections[3 * m2_idx + 0] -= phi_qm * lambda_spatial[0];
            result->force_corrections[3 * m2_idx + 1] -= phi_qm * lambda_spatial[1];
            result->force_corrections[3 * m2_idx + 2] -= phi_qm * lambda_spatial[2];

            /* Accumulate for M1 force */
            sum_phi_lambda[0] += phi_qm * lambda_spatial[0];
            sum_phi_lambda[1] += phi_qm * lambda_spatial[1];
            sum_phi_lambda[2] += phi_qm * lambda_spatial[2];
        }

        /* ThDD:T-US-035-6.10 — Force on M1: F_M1^(chain) = +sum_i Phi_QM(R_M2_i) * lambda_spatial */
        result->force_corrections[3 * m1_idx + 0] += sum_phi_lambda[0];
        result->force_corrections[3 * m1_idx + 1] += sum_phi_lambda[1];
        result->force_corrections[3 * m1_idx + 2] += sum_phi_lambda[2];
    } else {
        /* Full 4-constraint case: use lambda[1,2,3] directly */
        double lambda_spatial[3] = {lambda[1], lambda[2], lambda[3]};

        /* Sum of QM potential at M2 sites for M1 force */
        double sum_phi_lambda[3] = {0.0, 0.0, 0.0};

        for (int i = 0; i < n_neighbors; i++) {
            int m2_idx = m2_info->neighbors[i];
            const double *r_m2 = &charge_pos[3 * m2_idx];

            /* ThDD:T-US-035-6.5 — QM potential at M2 site */
            double phi_qm = compute_qm_potential(n_qm, qm_charges, qm_pos, r_m2);

            /* ThDD:T-US-035-6.9 — Force on M2: F_M2^(chain) = -Phi_QM * lambda */
            result->force_corrections[3 * m2_idx + 0] -= phi_qm * lambda_spatial[0];
            result->force_corrections[3 * m2_idx + 1] -= phi_qm * lambda_spatial[1];
            result->force_corrections[3 * m2_idx + 2] -= phi_qm * lambda_spatial[2];

            /* Accumulate for M1 force */
            sum_phi_lambda[0] += phi_qm * lambda_spatial[0];
            sum_phi_lambda[1] += phi_qm * lambda_spatial[1];
            sum_phi_lambda[2] += phi_qm * lambda_spatial[2];
        }

        /* ThDD:T-US-035-6.10 — Force on M1: F_M1^(chain) = +sum_i Phi_QM(R_M2_i) * lambda */
        result->force_corrections[3 * m1_idx + 0] += sum_phi_lambda[0];
        result->force_corrections[3 * m1_idx + 1] += sum_phi_lambda[1];
        result->force_corrections[3 * m1_idx + 2] += sum_phi_lambda[2];
    }

    /* Compute dipole error for diagnostics */
    /* ThDD:T-US-035-4.5 — sum(Delta Q * d) should be zero */
    double dipole_residual[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < n_neighbors; i++) {
        dipole_residual[0] += delta_q[i] * displacements[3 * i + 0];
        dipole_residual[1] += delta_q[i] * displacements[3 * i + 1];
        dipole_residual[2] += delta_q[i] * displacements[3 * i + 2];
    }
    result->dipole_error[0] += fabs(dipole_residual[0]);
    result->dipole_error[1] += fabs(dipole_residual[1]);
    result->dipole_error[2] += fabs(dipole_residual[2]);

    free(displacements);
    free(delta_q);

    return GRODFTB_SUCCESS;
}

/* ---------------------------------------------------------------------------
 * Internal: Apply shift scheme to all link atoms
 * ThDD:T-US-035-4.1-4.22 — Charge redistribution with dipole preservation
 * --------------------------------------------------------------------------- */
static int apply_shift_scheme(grodftb_linkatom_handle_t handle,
                              int ncharges,
                              const double *charges,
                              const double *charge_pos,
                              const double *qm_charges,
                              const double *qm_pos,
                              int n_qm,
                              grodftb_redistrib_result_t *result)
{
    /* Start with copy of original charges */
    memcpy(result->modified_charges, charges, ncharges * sizeof(double));

    /* Zero force corrections */
    memset(result->force_corrections, 0, 3 * ncharges * sizeof(double));

    /* Reset diagnostics */
    result->fallback_used = 0;
    result->charge_error = 0.0;
    result->dipole_error[0] = 0.0;
    result->dipole_error[1] = 0.0;
    result->dipole_error[2] = 0.0;

    /* Process each link atom */
    for (int i = 0; i < handle->nlinks; i++) {
        int link_fallback = 0;
        int rc = apply_shift_single_link(&handle->links[i],
                                         &handle->m2_info[i],
                                         ncharges,
                                         charge_pos,
                                         qm_pos,
                                         n_qm,
                                         qm_charges,
                                         result,
                                         &link_fallback);
        if (rc != GRODFTB_SUCCESS) {
            return rc;
        }
        if (link_fallback) {
            result->fallback_used = 1;
        }
    }

    /* Compute total charge conservation error */
    double sum_orig = 0.0;
    double sum_mod = 0.0;
    for (int i = 0; i < ncharges; i++) {
        sum_orig += charges[i];
        sum_mod += result->modified_charges[i];
    }
    result->charge_error = fabs(sum_mod - sum_orig);

    return GRODFTB_SUCCESS;
}

/* ---------------------------------------------------------------------------
 * Public: Main charge redistribution function
 * ThDD:T-US-035-2.1 (none), T-US-035-3.1 (zero), T-US-035-4.1-4.22 (shift)
 * --------------------------------------------------------------------------- */

int grodftb_linkatom_redistribute_charges(grodftb_linkatom_handle_t handle,
                                          int ncharges,
                                          const double *charges,
                                          const double *charge_pos,
                                          const double *qm_charges,
                                          const double *qm_pos,
                                          int n_qm,
                                          grodftb_redistrib_result_t *result)
{
    /* Validate inputs */
    if (handle == NULL) {
        return GRODFTB_ERR_NULL_HANDLE;
    }
    if (result == NULL) {
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    if (ncharges != result->ncharges) {
        return GRODFTB_ERR_SIZE_MISMATCH;
    }
    if (ncharges > 0 && (charges == NULL || charge_pos == NULL)) {
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    if (handle->charge_scheme == GRODFTB_CHARGE_SHIFT) {
        if (n_qm > 0 && (qm_charges == NULL || qm_pos == NULL)) {
            return GRODFTB_ERR_INVALID_ARGUMENT;
        }
    }

    /* Handle empty charge case */
    if (ncharges == 0) {
        result->fallback_used = 0;
        result->charge_error = 0.0;
        result->dipole_error[0] = 0.0;
        result->dipole_error[1] = 0.0;
        result->dipole_error[2] = 0.0;
        return GRODFTB_SUCCESS;
    }

    /* Dispatch to appropriate scheme */
    switch (handle->charge_scheme) {
    case GRODFTB_CHARGE_NONE:
        /* ThDD:T-US-035-2.1 — Passthrough */
        apply_none_scheme(ncharges, charges, result);
        break;

    case GRODFTB_CHARGE_ZERO:
        /* ThDD:T-US-035-3.1 — Zero M1 charges */
        apply_zero_scheme(ncharges, charges, handle->nlinks, handle->links, result);
        break;

    case GRODFTB_CHARGE_SHIFT:
        /* ThDD:T-US-035-4.1-4.22 — Full redistribution */
        return apply_shift_scheme(handle, ncharges, charges, charge_pos,
                                  qm_charges, qm_pos, n_qm, result);

    default:
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }

    return GRODFTB_SUCCESS;
}
