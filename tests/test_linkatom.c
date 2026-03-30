/*
 * SDD:specs.md:§9 — Unit tests for Module 5 (Link Atom Handler)
 * US-033: Link atom placement tests
 *
 * Tests validate:
 * - Link atom position placement (Eq. T-US-033-4.1)
 * - Link atom distance (Eq. T-US-033-2.3)
 * - Collinearity constraint (Eq. T-US-033-2.1)
 * - g-factor computation (Eq. T-US-033-4.4)
 * - Force conservation (Eq. T-US-033-5.4)
 * - Coordinate augmentation
 * - Species augmentation
 *
 * All tolerances from docs/theory/US-033/04_verification_criteria.md
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "grodftb/linkatom.h"
#include "grodftb/units.h"
#include "grodftb/error.h"

/* Test tolerances from 04_verification_criteria.md */
#define TOL_POSITION_NM     1e-10  /* Link position tolerance in nm */
#define TOL_POSITION_BOHR   (TOL_POSITION_NM * GRODFTB_NM_TO_BOHR)
#define TOL_COLLINEARITY    1e-12  /* Cross product magnitude */
#define TOL_FORCE_CONSERV   1e-12  /* Force conservation in atomic units */
#define TOL_G_FACTOR        1e-14  /* Relative error for g factor */

/* Test counters */
static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

#define RUN_TEST(test_func) do { \
    printf("  Running %s... ", #test_func); \
    tests_run++; \
    if (test_func()) { \
        printf("PASS\n"); \
        tests_passed++; \
    } else { \
        printf("FAIL\n"); \
        tests_failed++; \
    } \
} while(0)

/* ---------------------------------------------------------------------------
 * Helper functions
 * --------------------------------------------------------------------------- */

/* Compute Euclidean distance between two 3D vectors */
static double vec3_distance(const double *a, const double *b)
{
    double dx = b[0] - a[0];
    double dy = b[1] - a[1];
    double dz = b[2] - a[2];
    return sqrt(dx*dx + dy*dy + dz*dz);
}

/* Compute norm of 3D vector */
static double vec3_norm(const double *v)
{
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

/* Compute cross product c = a x b */
static void vec3_cross(const double *a, const double *b, double *c)
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

/* ---------------------------------------------------------------------------
 * Test: Link atom position for simple aligned geometry
 * ThDD:T-US-033-4.1 — R_L = R_A + d_link * unit_vec(R_B - R_A)
 * Config from docs/verification/US-033.md Appendix A.1
 * --------------------------------------------------------------------------- */
static int test_linkatom_position_simple(void)
{
    /*
     * Configuration: config_aligned_x
     * QM atom A: (0.0, 0.0, 0.0) Bohr
     * MM atom B: (2.910, 0.0, 0.0) Bohr  [0.154 nm = 2.910 Bohr, typical C-C]
     * d_link: 1.890 Bohr [0.1 nm]
     * Expected R_L: (1.890, 0.0, 0.0) Bohr
     */
    const double r_A[3] = {0.0, 0.0, 0.0};
    const double r_B[3] = {2.910, 0.0, 0.0};  /* 0.154 nm in Bohr */
    const double d_link = 1.890;              /* 0.1 nm in Bohr */
    const double expected_L[3] = {1.890, 0.0, 0.0};

    grodftb_linkatom_handle_t handle = NULL;
    int qm_atoms[1] = {0};
    int mm_atoms[1] = {0};  /* Global MM index, not used for position calc */
    double ref_lengths[1] = {d_link};
    int h_species = 0;

    int rc = grodftb_linkatom_create(1, qm_atoms, mm_atoms, ref_lengths,
                                     h_species, GRODFTB_CHARGE_NONE, &handle);
    if (rc != GRODFTB_SUCCESS) {
        printf("\n    ERROR: grodftb_linkatom_create failed with code %d\n", rc);
        return 0;
    }

    /* Create coordinate buffers for 1 QM atom and 1 MM atom */
    double qm_coords[3];
    memcpy(qm_coords, r_A, sizeof(qm_coords));

    double mm_coords[3];
    memcpy(mm_coords, r_B, sizeof(mm_coords));

    /* Compute link atom position */
    double link_pos[3] = {0.0, 0.0, 0.0};
    rc = grodftb_linkatom_compute_positions(handle, 1, qm_coords, mm_coords, link_pos);
    if (rc != GRODFTB_SUCCESS) {
        printf("\n    ERROR: grodftb_linkatom_compute_positions failed with code %d\n", rc);
        grodftb_linkatom_destroy(&handle);
        return 0;
    }

    /* Check position */
    double error = vec3_distance(link_pos, expected_L);
    if (error > TOL_POSITION_BOHR) {
        printf("\n    ERROR: Position error %.3e Bohr > tolerance %.3e\n",
               error, TOL_POSITION_BOHR);
        printf("    Expected: (%.6f, %.6f, %.6f)\n", expected_L[0], expected_L[1], expected_L[2]);
        printf("    Got:      (%.6f, %.6f, %.6f)\n", link_pos[0], link_pos[1], link_pos[2]);
        grodftb_linkatom_destroy(&handle);
        return 0;
    }

    grodftb_linkatom_destroy(&handle);
    return 1;
}

/* ---------------------------------------------------------------------------
 * Test: Link atom distance equals d_link exactly
 * ThDD:T-US-033-2.3 — d_AL = g * d_AB = d_link
 * --------------------------------------------------------------------------- */
static int test_linkatom_distance_fixed(void)
{
    /* Use arbitrary orientation */
    const double r_A[3] = {2.835, 3.779, 1.890};  /* (0.15, 0.20, 0.10) nm */
    const double r_B[3] = {5.669, 3.779, 1.890};  /* (0.30, 0.20, 0.10) nm */
    const double d_link = 1.890;                  /* 0.1 nm in Bohr */

    grodftb_linkatom_handle_t handle = NULL;
    int qm_atoms[1] = {0};
    int mm_atoms[1] = {0};
    double ref_lengths[1] = {d_link};
    int h_species = 0;

    int rc = grodftb_linkatom_create(1, qm_atoms, mm_atoms, ref_lengths,
                                     h_species, GRODFTB_CHARGE_NONE, &handle);
    if (rc != GRODFTB_SUCCESS) return 0;

    double qm_coords[3], mm_coords[3];
    memcpy(qm_coords, r_A, sizeof(qm_coords));
    memcpy(mm_coords, r_B, sizeof(mm_coords));

    double link_pos[3];
    rc = grodftb_linkatom_compute_positions(handle, 1, qm_coords, mm_coords, link_pos);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_linkatom_destroy(&handle);
        return 0;
    }

    /* Check d_AL = d_link */
    double d_AL = vec3_distance(link_pos, r_A);
    double error = fabs(d_AL - d_link);

    if (error > TOL_POSITION_BOHR) {
        printf("\n    ERROR: Distance error %.3e Bohr > tolerance %.3e\n",
               error, TOL_POSITION_BOHR);
        printf("    Expected d_AL: %.10f, Got: %.10f\n", d_link, d_AL);
        grodftb_linkatom_destroy(&handle);
        return 0;
    }

    grodftb_linkatom_destroy(&handle);
    return 1;
}

/* ---------------------------------------------------------------------------
 * Test: Link atom lies on A-B line (collinearity)
 * ThDD:T-US-033-2.1 — R_L = R_A + g*(R_B - R_A) implies L is on the A-B line
 * Verified by: (R_L - R_A) x (R_B - R_A) = 0
 * --------------------------------------------------------------------------- */
static int test_linkatom_collinearity(void)
{
    /* Use non-axis-aligned geometry */
    const double r_A[3] = {1.0, 2.0, 3.0};
    const double r_B[3] = {4.0, 5.0, 6.0};
    const double d_link = 1.5;

    grodftb_linkatom_handle_t handle = NULL;
    int qm_atoms[1] = {0};
    int mm_atoms[1] = {0};
    double ref_lengths[1] = {d_link};

    int rc = grodftb_linkatom_create(1, qm_atoms, mm_atoms, ref_lengths,
                                     0, GRODFTB_CHARGE_NONE, &handle);
    if (rc != GRODFTB_SUCCESS) return 0;

    double qm_coords[3], mm_coords[3];
    memcpy(qm_coords, r_A, sizeof(qm_coords));
    memcpy(mm_coords, r_B, sizeof(mm_coords));

    double link_pos[3];
    rc = grodftb_linkatom_compute_positions(handle, 1, qm_coords, mm_coords, link_pos);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_linkatom_destroy(&handle);
        return 0;
    }

    /* Compute (R_L - R_A) and (R_B - R_A) */
    double v1[3] = {link_pos[0] - r_A[0], link_pos[1] - r_A[1], link_pos[2] - r_A[2]};
    double v2[3] = {r_B[0] - r_A[0], r_B[1] - r_A[1], r_B[2] - r_A[2]};

    /* Cross product */
    double cross[3];
    vec3_cross(v1, v2, cross);
    double cross_norm = vec3_norm(cross);

    /* Normalize by product of norms for numerical stability */
    double norm_product = vec3_norm(v1) * vec3_norm(v2);
    double normalized_cross = (norm_product > 1e-30) ? cross_norm / norm_product : cross_norm;

    if (normalized_cross > TOL_COLLINEARITY) {
        printf("\n    ERROR: Collinearity violation %.3e > tolerance %.3e\n",
               normalized_cross, TOL_COLLINEARITY);
        grodftb_linkatom_destroy(&handle);
        return 0;
    }

    grodftb_linkatom_destroy(&handle);
    return 1;
}

/* ---------------------------------------------------------------------------
 * Test: g-factor computation
 * ThDD:T-US-033-4.4 — g = d_link / d_AB
 * --------------------------------------------------------------------------- */
static int test_linkatom_g_factor(void)
{
    const double r_A[3] = {0.0, 0.0, 0.0};
    const double r_B[3] = {3.0, 0.0, 0.0};  /* d_AB = 3.0 Bohr */
    const double d_link = 1.5;               /* d_link = 1.5 Bohr */
    const double expected_g = 0.5;           /* g = 1.5 / 3.0 = 0.5 */

    grodftb_linkatom_handle_t handle = NULL;
    int qm_atoms[1] = {0};
    int mm_atoms[1] = {0};
    double ref_lengths[1] = {d_link};

    int rc = grodftb_linkatom_create(1, qm_atoms, mm_atoms, ref_lengths,
                                     0, GRODFTB_CHARGE_NONE, &handle);
    if (rc != GRODFTB_SUCCESS) return 0;

    double qm_coords[3], mm_coords[3];
    memcpy(qm_coords, r_A, sizeof(qm_coords));
    memcpy(mm_coords, r_B, sizeof(mm_coords));

    /* Compute positions to update internal g factors */
    double link_pos[3];
    rc = grodftb_linkatom_compute_positions(handle, 1, qm_coords, mm_coords, link_pos);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_linkatom_destroy(&handle);
        return 0;
    }

    /* Get the g factor */
    double g_factor = 0.0;
    rc = grodftb_linkatom_get_g_factor(handle, 0, &g_factor);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_linkatom_destroy(&handle);
        return 0;
    }

    double rel_error = fabs(g_factor - expected_g) / expected_g;
    if (rel_error > TOL_G_FACTOR) {
        printf("\n    ERROR: g-factor error %.3e > tolerance %.3e\n",
               rel_error, TOL_G_FACTOR);
        printf("    Expected: %.15f, Got: %.15f\n", expected_g, g_factor);
        grodftb_linkatom_destroy(&handle);
        return 0;
    }

    grodftb_linkatom_destroy(&handle);
    return 1;
}

/* ---------------------------------------------------------------------------
 * Test: Force conservation F_A + F_B = F_L
 * ThDD:T-US-033-5.4 — Force projection must conserve total force
 * --------------------------------------------------------------------------- */
static int test_linkatom_force_conservation(void)
{
    /* Setup geometry */
    const double r_A[3] = {0.0, 0.0, 0.0};
    const double r_B[3] = {3.0, 0.0, 0.0};
    const double d_link = 1.5;

    /* Arbitrary force on link atom */
    const double F_L[3] = {0.01, 0.02, -0.015};

    grodftb_linkatom_handle_t handle = NULL;
    int qm_atoms[1] = {0};
    int mm_atoms[1] = {0};
    double ref_lengths[1] = {d_link};

    int rc = grodftb_linkatom_create(1, qm_atoms, mm_atoms, ref_lengths,
                                     0, GRODFTB_CHARGE_NONE, &handle);
    if (rc != GRODFTB_SUCCESS) return 0;

    double qm_coords[3], mm_coords[3];
    memcpy(qm_coords, r_A, sizeof(qm_coords));
    memcpy(mm_coords, r_B, sizeof(mm_coords));

    /* First compute positions to set up g factors */
    double link_pos[3];
    rc = grodftb_linkatom_compute_positions(handle, 1, qm_coords, mm_coords, link_pos);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_linkatom_destroy(&handle);
        return 0;
    }

    /* Create augmented force array: [F_QM0, F_L0] */
    double augmented_forces[6] = {0.0, 0.0, 0.0,  /* F on QM atom (initially 0) */
                                  F_L[0], F_L[1], F_L[2]};  /* F on link atom */

    /* Output force arrays */
    double qm_forces[3] = {0.0, 0.0, 0.0};
    double mm_forces[3] = {0.0, 0.0, 0.0};

    rc = grodftb_linkatom_project_forces(handle, 1, qm_coords, mm_coords,
                                         augmented_forces, qm_forces, mm_forces);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_linkatom_destroy(&handle);
        return 0;
    }

    /* Check F_A + F_B = F_L */
    double F_sum[3] = {qm_forces[0] + mm_forces[0],
                       qm_forces[1] + mm_forces[1],
                       qm_forces[2] + mm_forces[2]};

    double error[3] = {F_sum[0] - F_L[0], F_sum[1] - F_L[1], F_sum[2] - F_L[2]};
    double error_norm = vec3_norm(error);

    if (error_norm > TOL_FORCE_CONSERV) {
        printf("\n    ERROR: Force conservation violation %.3e > tolerance %.3e\n",
               error_norm, TOL_FORCE_CONSERV);
        printf("    F_A: (%.6e, %.6e, %.6e)\n", qm_forces[0], qm_forces[1], qm_forces[2]);
        printf("    F_B: (%.6e, %.6e, %.6e)\n", mm_forces[0], mm_forces[1], mm_forces[2]);
        printf("    Sum: (%.6e, %.6e, %.6e)\n", F_sum[0], F_sum[1], F_sum[2]);
        printf("    F_L: (%.6e, %.6e, %.6e)\n", F_L[0], F_L[1], F_L[2]);
        grodftb_linkatom_destroy(&handle);
        return 0;
    }

    grodftb_linkatom_destroy(&handle);
    return 1;
}

/* ---------------------------------------------------------------------------
 * Test: Coordinate array augmentation
 * SDD:specs.md:§9.5 — Augmented array layout: [QM_0, ..., QM_{n-1}, L_0, ...]
 * --------------------------------------------------------------------------- */
static int test_linkatom_augment_coords(void)
{
    /* 2 QM atoms, 1 link atom */
    const double qm_coords[6] = {0.0, 0.0, 0.0,   /* QM atom 0 */
                                 1.0, 2.0, 3.0};  /* QM atom 1 */
    const double mm_coords[3] = {3.0, 0.0, 0.0};  /* MM atom (boundary partner of QM atom 0) */
    const double d_link = 1.5;

    grodftb_linkatom_handle_t handle = NULL;
    int qm_atoms[1] = {0};  /* Link from QM atom 0 */
    int mm_atoms[1] = {0};
    double ref_lengths[1] = {d_link};

    int rc = grodftb_linkatom_create(1, qm_atoms, mm_atoms, ref_lengths,
                                     0, GRODFTB_CHARGE_NONE, &handle);
    if (rc != GRODFTB_SUCCESS) return 0;

    /* Augmented array: 2 QM atoms + 1 link atom = 9 doubles */
    double augmented[9];
    rc = grodftb_linkatom_augment_coords(handle, 2, qm_coords, mm_coords, augmented);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_linkatom_destroy(&handle);
        return 0;
    }

    /* Check QM atoms preserved */
    int ok = 1;
    for (int i = 0; i < 6; i++) {
        if (fabs(augmented[i] - qm_coords[i]) > 1e-15) {
            printf("\n    ERROR: QM coord[%d] not preserved: %.15f != %.15f\n",
                   i, augmented[i], qm_coords[i]);
            ok = 0;
        }
    }

    /* Check link atom at expected position (indices 6,7,8) */
    /* Expected link: (0,0,0) + (1.5/3.0) * (3,0,0) - (0,0,0) = (1.5, 0, 0) */
    const double expected_link[3] = {1.5, 0.0, 0.0};
    double link_error = vec3_distance(&augmented[6], expected_link);
    if (link_error > TOL_POSITION_BOHR) {
        printf("\n    ERROR: Link atom position error %.3e > tolerance\n", link_error);
        printf("    Expected: (%.6f, %.6f, %.6f)\n", expected_link[0], expected_link[1], expected_link[2]);
        printf("    Got:      (%.6f, %.6f, %.6f)\n", augmented[6], augmented[7], augmented[8]);
        ok = 0;
    }

    grodftb_linkatom_destroy(&handle);
    return ok;
}

/* ---------------------------------------------------------------------------
 * Test: Species array augmentation
 * SDD:specs.md:§9.2 — Link atoms assigned H species
 * --------------------------------------------------------------------------- */
static int test_linkatom_species_augmentation(void)
{
    /* 2 QM atoms with species C and N, 1 link atom */
    const int qm_species[2] = {1, 2};  /* C=1, N=2 */
    const int h_species = 0;           /* H=0 */

    grodftb_linkatom_handle_t handle = NULL;
    int qm_atoms[1] = {0};
    int mm_atoms[1] = {0};
    double ref_lengths[1] = {1.5};

    int rc = grodftb_linkatom_create(1, qm_atoms, mm_atoms, ref_lengths,
                                     h_species, GRODFTB_CHARGE_NONE, &handle);
    if (rc != GRODFTB_SUCCESS) return 0;

    /* Augmented species: 2 QM atoms + 1 link atom = 3 integers */
    int augmented_species[3];
    rc = grodftb_linkatom_augment_species(handle, 2, qm_species, augmented_species);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_linkatom_destroy(&handle);
        return 0;
    }

    /* Check QM species preserved */
    int ok = 1;
    if (augmented_species[0] != qm_species[0]) {
        printf("\n    ERROR: QM species[0] not preserved: %d != %d\n",
               augmented_species[0], qm_species[0]);
        ok = 0;
    }
    if (augmented_species[1] != qm_species[1]) {
        printf("\n    ERROR: QM species[1] not preserved: %d != %d\n",
               augmented_species[1], qm_species[1]);
        ok = 0;
    }

    /* Check link atom has H species */
    if (augmented_species[2] != h_species) {
        printf("\n    ERROR: Link atom species wrong: %d != %d (H)\n",
               augmented_species[2], h_species);
        ok = 0;
    }

    grodftb_linkatom_destroy(&handle);
    return ok;
}

/* ---------------------------------------------------------------------------
 * Test: Error on zero distance (collapsed atoms)
 * --------------------------------------------------------------------------- */
static int test_linkatom_zero_distance_error(void)
{
    const double r_A[3] = {1.0, 2.0, 3.0};
    const double r_B[3] = {1.0, 2.0, 3.0};  /* Same as A! */
    const double d_link = 1.5;

    grodftb_linkatom_handle_t handle = NULL;
    int qm_atoms[1] = {0};
    int mm_atoms[1] = {0};
    double ref_lengths[1] = {d_link};

    int rc = grodftb_linkatom_create(1, qm_atoms, mm_atoms, ref_lengths,
                                     0, GRODFTB_CHARGE_NONE, &handle);
    if (rc != GRODFTB_SUCCESS) return 0;

    double qm_coords[3], mm_coords[3];
    memcpy(qm_coords, r_A, sizeof(qm_coords));
    memcpy(mm_coords, r_B, sizeof(mm_coords));

    double link_pos[3];
    rc = grodftb_linkatom_compute_positions(handle, 1, qm_coords, mm_coords, link_pos);

    /* Should return an error */
    if (rc == GRODFTB_SUCCESS) {
        printf("\n    ERROR: Expected error for zero distance, got success\n");
        grodftb_linkatom_destroy(&handle);
        return 0;
    }

    grodftb_linkatom_destroy(&handle);
    return 1;  /* Error correctly detected */
}

/* ---------------------------------------------------------------------------
 * Test: Get link count
 * --------------------------------------------------------------------------- */
static int test_linkatom_count(void)
{
    grodftb_linkatom_handle_t handle = NULL;
    int qm_atoms[3] = {0, 1, 2};
    int mm_atoms[3] = {5, 6, 7};
    double ref_lengths[3] = {1.5, 1.5, 1.5};

    int rc = grodftb_linkatom_create(3, qm_atoms, mm_atoms, ref_lengths,
                                     0, GRODFTB_CHARGE_NONE, &handle);
    if (rc != GRODFTB_SUCCESS) return 0;

    int count = grodftb_linkatom_count(handle);
    if (count != 3) {
        printf("\n    ERROR: Expected count 3, got %d\n", count);
        grodftb_linkatom_destroy(&handle);
        return 0;
    }

    grodftb_linkatom_destroy(&handle);
    return 1;
}

/* ---------------------------------------------------------------------------
 * Test: Null pointer handling
 * --------------------------------------------------------------------------- */
static int test_linkatom_null_pointers(void)
{
    /* NULL handle */
    int count = grodftb_linkatom_count(NULL);
    if (count != 0) {
        printf("\n    ERROR: Expected 0 for NULL handle, got %d\n", count);
        return 0;
    }

    /* NULL output in create */
    int rc = grodftb_linkatom_create(0, NULL, NULL, NULL, 0, 0, NULL);
    if (rc == GRODFTB_SUCCESS) {
        printf("\n    ERROR: Expected error for NULL output, got success\n");
        return 0;
    }

    return 1;
}

/* ---------------------------------------------------------------------------
 * Main test runner
 * --------------------------------------------------------------------------- */
int main(void)
{
    printf("=== US-033 Link Atom Unit Tests ===\n\n");

    printf("Geometric verification:\n");
    RUN_TEST(test_linkatom_position_simple);
    RUN_TEST(test_linkatom_distance_fixed);
    RUN_TEST(test_linkatom_collinearity);
    RUN_TEST(test_linkatom_g_factor);

    printf("\nForce projection:\n");
    RUN_TEST(test_linkatom_force_conservation);

    printf("\nArray augmentation:\n");
    RUN_TEST(test_linkatom_augment_coords);
    RUN_TEST(test_linkatom_species_augmentation);

    printf("\nEdge cases:\n");
    RUN_TEST(test_linkatom_zero_distance_error);
    RUN_TEST(test_linkatom_count);
    RUN_TEST(test_linkatom_null_pointers);

    printf("\n=== Summary ===\n");
    printf("Tests run:    %d\n", tests_run);
    printf("Tests passed: %d\n", tests_passed);
    printf("Tests failed: %d\n", tests_failed);

    return tests_failed > 0 ? 1 : 0;
}
