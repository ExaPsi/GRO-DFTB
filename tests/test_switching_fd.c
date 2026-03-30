/*
 * SDD:specs.md:S21.1 -- Finite-difference force validation for switched embedding
 * US-039: Switching Function for Cutoff Embedding
 *
 * ThDD:T-US-039-V.4 -- FD validation for MM forces in switching region
 *
 * These tests validate that analytic forces on MM atoms are consistent with
 * the energy gradient computed by finite differences when the switching
 * function is active.
 *
 * IMPORTANT: Both the analytic forces and the FD forces are computed during
 * the test. Neither side uses assumed or fabricated reference values.
 * This is the gold standard for force validation per CLAUDE.md Golden Rule.
 *
 * FD methodology (ThDD:T-US-038-V.2):
 *   - Central differences: F = -(E(+d) - E(-d)) / (2*d)
 *   - Step size: delta = 10^-4 Bohr
 *   - Tolerance: 10^-4 relative error (specs.md S21.1)
 *   - Absolute floor: 10^-6 Ha/Bohr (for small force components)
 *
 * LGPL-3.0-or-later
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "grodftb/driver.h"
#include "grodftb/error.h"
#include "grodftb/units.h"

/* ---------------------------------------------------------------------------
 * FD Test Parameters (from docs/theory/US-038/03_numerical_methods.md)
 *
 * ThDD:T-US-038-V.2 -- Recommended FD parameters, also applicable to US-039
 * --------------------------------------------------------------------------- */
#define FD_DELTA            1e-4   /* Bohr -- ThDD:T-US-038-V.2 */
#define FD_REL_TOLERANCE    1e-4   /* Relative error tolerance */
#define FD_ABS_FLOOR        1e-6   /* Ha/Bohr -- absolute tolerance floor */

/* ---------------------------------------------------------------------------
 * Test framework
 * --------------------------------------------------------------------------- */
static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

#define RUN_TEST(test_func) do { \
    printf("  Running %s... ", #test_func); \
    fflush(stdout); \
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
 * B4 benchmark geometry (from tests/data/b4/)
 *
 * Provenance: tests/data/b4/provenance.txt
 *   - DFTB+ 25.1 (commit fd31d873)
 *   - SK parameters: mio-1-1
 * --------------------------------------------------------------------------- */

/* B1/B4 QM geometry (water dimer) in Bohr
 * Source: tests/data/b1/geo.gen */
static const double B1_COORDS_BOHR[18] = {
    -1.326958030291, -0.105938038921,  0.018787655779,  /* O1 */
    -1.931664677465,  1.600174613723, -0.021711061883,  /* H1 */
     0.486644126310,  0.079597148366,  0.009862479935,  /* H2 */
     4.196837646028,  0.050487809237,  0.001171630113,  /* O2 */
     4.908550027307, -0.777930269645,  1.448937953129,  /* H3 */
     4.900314601448, -0.449424272972, -1.407433901241   /* H4 */
};

/* Species indices: O=0, H=1 */
static const int B1_SPECIES[6] = { 0, 1, 1, 0, 1, 1 };

#define B4_N_QM    6
#define B4_N_MM    1
#define B4_MM_CHARGE  1.0

/* Conversion factor: Angstrom to Bohr (CODATA 2022) */
#define ANGSTROM_TO_BOHR  1.8897259886

/* ---------------------------------------------------------------------------
 * Helper: Path resolution for test data
 * --------------------------------------------------------------------------- */
static const char *get_b1_hsd_path(void)
{
    static char path[1024];
    const char *srcdir = getenv("GRODFTB_SOURCE_DIR");
    if (!srcdir) {
        srcdir = GRODFTB_SOURCE_DIR_STR;
    }
    snprintf(path, sizeof(path), "%s/tests/data/b1/dftb_in.hsd", srcdir);
    return path;
}

static const char *get_b1_dir(void)
{
    static char path[1024];
    const char *srcdir = getenv("GRODFTB_SOURCE_DIR");
    if (!srcdir) {
        srcdir = GRODFTB_SOURCE_DIR_STR;
    }
    snprintf(path, sizeof(path), "%s/tests/data/b1", srcdir);
    return path;
}

/* ---------------------------------------------------------------------------
 * ThDD:T-US-038-V.2 -- Combined acceptance criterion
 *
 * The test passes if:
 *   |F_analytic - F_FD| < max(rel_tol * |F_analytic|, abs_floor)
 * --------------------------------------------------------------------------- */
static int check_fd_tolerance(double f_analytic, double f_fd,
                               double *abs_error_out, double *rel_error_out)
{
    double abs_error = fabs(f_analytic - f_fd);
    double threshold = fmax(FD_REL_TOLERANCE * fabs(f_analytic), FD_ABS_FLOOR);

    if (abs_error_out) *abs_error_out = abs_error;
    if (rel_error_out) {
        *rel_error_out = (fabs(f_analytic) > 1e-15)
                         ? abs_error / fabs(f_analytic)
                         : abs_error;
    }

    return (abs_error < threshold) ? 1 : 0;
}

/* ---------------------------------------------------------------------------
 * Test 1: test_fd_mm_forces_switching_mid_T_US_039_V_4
 *
 * ThDD:T-US-039-V.4 -- FD validation at u = 0.5 (middle of switching region)
 *
 * This tests the point where S'(u) is maximal in magnitude:
 *   S'(0.5) = -30 * (0.5)^2 * (0.5)^2 = -30/16 = -1.875
 *
 * The switching derivative contributes maximally to the force at this point,
 * so it's the most demanding test for the switched force formula.
 *
 * Test setup:
 *   - QM region: B4 water dimer
 *   - MM charge: +1.0 e at distance r such that u = 0.5
 *   - Switching: r_on = 5.0 Bohr, switch_width = 2.0 Bohr, r_off = 7.0 Bohr
 *   - MM position: r = r_on + 0.5 * switch_width = 6.0 Bohr from QM center
 *
 * IMPORTANT: Both analytic and FD forces are computed, not assumed.
 * --------------------------------------------------------------------------- */
static int test_fd_mm_forces_switching_mid_T_US_039_V_4(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    printf(" SKIPPED (no DFTB+)");
    return 1;
#else
    const char *datadir = get_b1_dir();
    const char *hsd_path = get_b1_hsd_path();
    char oldcwd[1024];

    if (!getcwd(oldcwd, sizeof(oldcwd))) {
        printf("\n    getcwd failed\n");
        return 0;
    }
    if (chdir(datadir) != 0) {
        printf("\n    chdir to b1 data dir failed\n");
        return 0;
    }

    /*
     * Switching region parameters (in Bohr):
     *   r_on = 5.0 Bohr (switching starts)
     *   switch_width = 2.0 Bohr
     *   r_off = 7.0 Bohr (cutoff)
     *
     * For u = 0.5: r = r_on + 0.5 * switch_width = 6.0 Bohr
     *
     * Place MM charge at (-6.0, 0.0, 0.0) Bohr from the QM center of mass.
     * QM center of mass is approximately at (1.2, 0.0, 0.0) Bohr.
     * So absolute position: (-6.0 + 1.2, 0.0, 0.0) = (-4.8, 0.0, 0.0) Bohr.
     *
     * For simplicity, use a position that puts the charge at ~ 6 Bohr from O1.
     */
    double mm_pos[3] = { -7.0, 0.0, 0.0 };  /* ~6 Bohr from O1 at (-1.3, 0, 0) */
    double mm_charge = B4_MM_CHARGE;

    /* Switching parameters to pass to embedding API (TO BE IMPLEMENTED) */
    double switch_width_bohr = 2.0;
    double cutoff_bohr = 7.0;

    /*
     * First, compute analytic forces with switching enabled
     */
    grodftb_handle_t handle = NULL;
    int rc = grodftb_init(hsd_path, &handle);
    if (rc != GRODFTB_SUCCESS) {
        chdir(oldcwd);
        printf("\n    grodftb_init failed: %d\n", rc);
        return 0;
    }

    rc = grodftb_set_geometry(handle, B4_N_QM, B1_SPECIES, B1_COORDS_BOHR);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        chdir(oldcwd);
        printf("\n    grodftb_set_geometry failed: %d\n", rc);
        return 0;
    }

    /*
     * Set embedding with switching parameters (API TO BE IMPLEMENTED)
     *
     * Expected API:
     *   grodftb_set_embedding_params(handle, cutoff_bohr, switch_width_bohr);
     *   grodftb_set_embedding_charges(handle, n_mm, charges, positions);
     *
     * For now, this will fail because the API doesn't exist yet (TDD RED phase).
     */
    rc = grodftb_set_embedding_params(handle, cutoff_bohr, switch_width_bohr);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        chdir(oldcwd);
        printf("\n    grodftb_set_embedding_params failed: %d (expected in RED phase)\n", rc);
        return 0;
    }

    rc = grodftb_set_embedding_charges(handle, B4_N_MM, &mm_charge, mm_pos);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        chdir(oldcwd);
        printf("\n    grodftb_set_embedding_charges failed: %d\n", rc);
        return 0;
    }

    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        chdir(oldcwd);
        printf("\n    grodftb_compute failed: %d\n", rc);
        return 0;
    }

    /* Get analytic forces on MM charges */
    double f_analytic[3];
    rc = grodftb_get_embedding_forces(handle, f_analytic);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        chdir(oldcwd);
        printf("\n    grodftb_get_embedding_forces failed: %d\n", rc);
        return 0;
    }

    grodftb_finalize(&handle);

    /*
     * Compute FD forces by perturbing MM position
     *
     * ThDD:T-US-039-V.4 -- Central difference formula:
     *   F_alpha = -(E(R_J + delta*e_alpha) - E(R_J - delta*e_alpha)) / (2*delta)
     */
    double f_fd[3];
    int all_pass = 1;
    const char *dir_names[] = {"x", "y", "z"};

    printf("\n");
    printf("    Switching region: r_on=%.1f Bohr, width=%.1f Bohr, r_off=%.1f Bohr\n",
           cutoff_bohr - switch_width_bohr, switch_width_bohr, cutoff_bohr);
    printf("    MM position: (%.1f, %.1f, %.1f) Bohr\n",
           mm_pos[0], mm_pos[1], mm_pos[2]);
    printf("    Testing at u ~ 0.5 (maximum S' contribution)\n\n");

    for (int dir = 0; dir < 3; dir++) {
        double e_plus = 0.0, e_minus = 0.0;

        /* Positive displacement */
        {
            double mm_pos_pert[3] = {mm_pos[0], mm_pos[1], mm_pos[2]};
            mm_pos_pert[dir] += FD_DELTA;

            handle = NULL;
            rc = grodftb_init(hsd_path, &handle);
            if (rc != GRODFTB_SUCCESS) {
                chdir(oldcwd);
                printf("    grodftb_init failed for +delta\n");
                return 0;
            }

            grodftb_set_geometry(handle, B4_N_QM, B1_SPECIES, B1_COORDS_BOHR);
            grodftb_set_embedding_params(handle, cutoff_bohr, switch_width_bohr);
            grodftb_set_embedding_charges(handle, B4_N_MM, &mm_charge, mm_pos_pert);
            grodftb_compute(handle);
            grodftb_get_energy(handle, &e_plus);
            grodftb_finalize(&handle);
        }

        /* Negative displacement */
        {
            double mm_pos_pert[3] = {mm_pos[0], mm_pos[1], mm_pos[2]};
            mm_pos_pert[dir] -= FD_DELTA;

            handle = NULL;
            rc = grodftb_init(hsd_path, &handle);
            if (rc != GRODFTB_SUCCESS) {
                chdir(oldcwd);
                printf("    grodftb_init failed for -delta\n");
                return 0;
            }

            grodftb_set_geometry(handle, B4_N_QM, B1_SPECIES, B1_COORDS_BOHR);
            grodftb_set_embedding_params(handle, cutoff_bohr, switch_width_bohr);
            grodftb_set_embedding_charges(handle, B4_N_MM, &mm_charge, mm_pos_pert);
            grodftb_compute(handle);
            grodftb_get_energy(handle, &e_minus);
            grodftb_finalize(&handle);
        }

        /* ThDD:T-US-038-V.2 -- Central difference */
        f_fd[dir] = -(e_plus - e_minus) / (2.0 * FD_DELTA);

        /* Check tolerance */
        double abs_err, rel_err;
        int pass = check_fd_tolerance(f_analytic[dir], f_fd[dir], &abs_err, &rel_err);

        printf("    F_MM[%s]: analytic=%.10e  FD=%.10e  rel_err=%.2e  %s\n",
               dir_names[dir], f_analytic[dir], f_fd[dir], rel_err,
               pass ? "OK" : "FAIL");

        if (!pass) all_pass = 0;
    }

    chdir(oldcwd);
    return all_pass;
#endif
}

/* ---------------------------------------------------------------------------
 * Test 2: test_fd_mm_forces_switching_early_T_US_039_V_4
 *
 * ThDD:T-US-039-V.4 -- FD validation at u = 0.25 (early switching region)
 *
 * At this point:
 *   S(0.25) ~ 0.896
 *   S'(0.25) = -30 * (0.25)^2 * (0.75)^2 = -1.0546875
 *
 * The force is dominated by the unswitched Coulomb term with small S' correction.
 * --------------------------------------------------------------------------- */
static int test_fd_mm_forces_switching_early_T_US_039_V_4(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    printf(" SKIPPED (no DFTB+)");
    return 1;
#else
    const char *datadir = get_b1_dir();
    const char *hsd_path = get_b1_hsd_path();
    char oldcwd[1024];

    if (!getcwd(oldcwd, sizeof(oldcwd))) {
        printf("\n    getcwd failed\n");
        return 0;
    }
    if (chdir(datadir) != 0) {
        printf("\n    chdir to b1 data dir failed\n");
        return 0;
    }

    /* For u = 0.25: r = r_on + 0.25 * switch_width = 5.0 + 0.5 = 5.5 Bohr */
    double mm_pos[3] = { -6.5, 0.0, 0.0 };  /* ~ 5.5 Bohr from O1 */
    double mm_charge = B4_MM_CHARGE;

    double switch_width_bohr = 2.0;
    double cutoff_bohr = 7.0;

    /* Get analytic forces */
    grodftb_handle_t handle = NULL;
    int rc = grodftb_init(hsd_path, &handle);
    if (rc != GRODFTB_SUCCESS) {
        chdir(oldcwd);
        printf("\n    grodftb_init failed: %d\n", rc);
        return 0;
    }

    grodftb_set_geometry(handle, B4_N_QM, B1_SPECIES, B1_COORDS_BOHR);
    rc = grodftb_set_embedding_params(handle, cutoff_bohr, switch_width_bohr);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        chdir(oldcwd);
        printf("\n    grodftb_set_embedding_params failed: %d\n", rc);
        return 0;
    }
    grodftb_set_embedding_charges(handle, B4_N_MM, &mm_charge, mm_pos);
    grodftb_compute(handle);

    double f_analytic[3];
    rc = grodftb_get_embedding_forces(handle, f_analytic);
    grodftb_finalize(&handle);

    if (rc != GRODFTB_SUCCESS) {
        chdir(oldcwd);
        printf("\n    grodftb_get_embedding_forces failed: %d\n", rc);
        return 0;
    }

    /* Compute FD forces */
    double f_fd[3];
    int all_pass = 1;
    const char *dir_names[] = {"x", "y", "z"};

    printf("\n");
    printf("    Testing at u ~ 0.25 (early switching)\n\n");

    for (int dir = 0; dir < 3; dir++) {
        double e_plus = 0.0, e_minus = 0.0;

        /* +delta */
        {
            double mm_pos_pert[3] = {mm_pos[0], mm_pos[1], mm_pos[2]};
            mm_pos_pert[dir] += FD_DELTA;

            handle = NULL;
            grodftb_init(hsd_path, &handle);
            grodftb_set_geometry(handle, B4_N_QM, B1_SPECIES, B1_COORDS_BOHR);
            grodftb_set_embedding_params(handle, cutoff_bohr, switch_width_bohr);
            grodftb_set_embedding_charges(handle, B4_N_MM, &mm_charge, mm_pos_pert);
            grodftb_compute(handle);
            grodftb_get_energy(handle, &e_plus);
            grodftb_finalize(&handle);
        }

        /* -delta */
        {
            double mm_pos_pert[3] = {mm_pos[0], mm_pos[1], mm_pos[2]};
            mm_pos_pert[dir] -= FD_DELTA;

            handle = NULL;
            grodftb_init(hsd_path, &handle);
            grodftb_set_geometry(handle, B4_N_QM, B1_SPECIES, B1_COORDS_BOHR);
            grodftb_set_embedding_params(handle, cutoff_bohr, switch_width_bohr);
            grodftb_set_embedding_charges(handle, B4_N_MM, &mm_charge, mm_pos_pert);
            grodftb_compute(handle);
            grodftb_get_energy(handle, &e_minus);
            grodftb_finalize(&handle);
        }

        f_fd[dir] = -(e_plus - e_minus) / (2.0 * FD_DELTA);

        double abs_err, rel_err;
        int pass = check_fd_tolerance(f_analytic[dir], f_fd[dir], &abs_err, &rel_err);

        printf("    F_MM[%s]: analytic=%.10e  FD=%.10e  rel_err=%.2e  %s\n",
               dir_names[dir], f_analytic[dir], f_fd[dir], rel_err,
               pass ? "OK" : "FAIL");

        if (!pass) all_pass = 0;
    }

    chdir(oldcwd);
    return all_pass;
#endif
}

/* ---------------------------------------------------------------------------
 * Test 3: test_fd_mm_forces_switching_late_T_US_039_V_4
 *
 * ThDD:T-US-039-V.4 -- FD validation at u = 0.75 (late switching region)
 *
 * At this point:
 *   S(0.75) ~ 0.104
 *   S'(0.75) = -30 * (0.75)^2 * (0.25)^2 = -1.0546875
 *
 * The force is nearly zero (small S) with significant S' correction.
 * This tests the force formula near the cutoff boundary.
 * --------------------------------------------------------------------------- */
static int test_fd_mm_forces_switching_late_T_US_039_V_4(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    printf(" SKIPPED (no DFTB+)");
    return 1;
#else
    const char *datadir = get_b1_dir();
    const char *hsd_path = get_b1_hsd_path();
    char oldcwd[1024];

    if (!getcwd(oldcwd, sizeof(oldcwd))) {
        printf("\n    getcwd failed\n");
        return 0;
    }
    if (chdir(datadir) != 0) {
        printf("\n    chdir to b1 data dir failed\n");
        return 0;
    }

    /* For u = 0.75: r = r_on + 0.75 * switch_width = 5.0 + 1.5 = 6.5 Bohr */
    double mm_pos[3] = { -7.5, 0.0, 0.0 };  /* ~ 6.5 Bohr from O1 */
    double mm_charge = B4_MM_CHARGE;

    double switch_width_bohr = 2.0;
    double cutoff_bohr = 7.0;

    /* Get analytic forces */
    grodftb_handle_t handle = NULL;
    int rc = grodftb_init(hsd_path, &handle);
    if (rc != GRODFTB_SUCCESS) {
        chdir(oldcwd);
        printf("\n    grodftb_init failed: %d\n", rc);
        return 0;
    }

    grodftb_set_geometry(handle, B4_N_QM, B1_SPECIES, B1_COORDS_BOHR);
    rc = grodftb_set_embedding_params(handle, cutoff_bohr, switch_width_bohr);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        chdir(oldcwd);
        printf("\n    grodftb_set_embedding_params failed: %d\n", rc);
        return 0;
    }
    grodftb_set_embedding_charges(handle, B4_N_MM, &mm_charge, mm_pos);
    grodftb_compute(handle);

    double f_analytic[3];
    rc = grodftb_get_embedding_forces(handle, f_analytic);
    grodftb_finalize(&handle);

    if (rc != GRODFTB_SUCCESS) {
        chdir(oldcwd);
        printf("\n    grodftb_get_embedding_forces failed: %d\n", rc);
        return 0;
    }

    /* Compute FD forces */
    double f_fd[3];
    int all_pass = 1;
    const char *dir_names[] = {"x", "y", "z"};

    printf("\n");
    printf("    Testing at u ~ 0.75 (late switching)\n\n");

    for (int dir = 0; dir < 3; dir++) {
        double e_plus = 0.0, e_minus = 0.0;

        /* +delta */
        {
            double mm_pos_pert[3] = {mm_pos[0], mm_pos[1], mm_pos[2]};
            mm_pos_pert[dir] += FD_DELTA;

            handle = NULL;
            grodftb_init(hsd_path, &handle);
            grodftb_set_geometry(handle, B4_N_QM, B1_SPECIES, B1_COORDS_BOHR);
            grodftb_set_embedding_params(handle, cutoff_bohr, switch_width_bohr);
            grodftb_set_embedding_charges(handle, B4_N_MM, &mm_charge, mm_pos_pert);
            grodftb_compute(handle);
            grodftb_get_energy(handle, &e_plus);
            grodftb_finalize(&handle);
        }

        /* -delta */
        {
            double mm_pos_pert[3] = {mm_pos[0], mm_pos[1], mm_pos[2]};
            mm_pos_pert[dir] -= FD_DELTA;

            handle = NULL;
            grodftb_init(hsd_path, &handle);
            grodftb_set_geometry(handle, B4_N_QM, B1_SPECIES, B1_COORDS_BOHR);
            grodftb_set_embedding_params(handle, cutoff_bohr, switch_width_bohr);
            grodftb_set_embedding_charges(handle, B4_N_MM, &mm_charge, mm_pos_pert);
            grodftb_compute(handle);
            grodftb_get_energy(handle, &e_minus);
            grodftb_finalize(&handle);
        }

        f_fd[dir] = -(e_plus - e_minus) / (2.0 * FD_DELTA);

        double abs_err, rel_err;
        int pass = check_fd_tolerance(f_analytic[dir], f_fd[dir], &abs_err, &rel_err);

        printf("    F_MM[%s]: analytic=%.10e  FD=%.10e  rel_err=%.2e  %s\n",
               dir_names[dir], f_analytic[dir], f_fd[dir], rel_err,
               pass ? "OK" : "FAIL");

        if (!pass) all_pass = 0;
    }

    chdir(oldcwd);
    return all_pass;
#endif
}

/* ---------------------------------------------------------------------------
 * main
 * --------------------------------------------------------------------------- */
int main(void)
{
    printf("=== US-039 Switching FD Force Validation Tests ===\n\n");
    printf("ThDD:T-US-039-V.4 -- FD validation for MM forces in switching region\n");
    printf("SDD:specs.md:S21.1 -- Tolerance: 10^-4 relative error\n\n");

    printf("Test methodology:\n");
    printf("  - Central differences: F = -(E(+d) - E(-d)) / (2*d)\n");
    printf("  - Step size: delta = 10^-4 Bohr\n");
    printf("  - Both analytic and FD forces are computed (no fabricated values)\n\n");

    RUN_TEST(test_fd_mm_forces_switching_mid_T_US_039_V_4);
    RUN_TEST(test_fd_mm_forces_switching_early_T_US_039_V_4);
    RUN_TEST(test_fd_mm_forces_switching_late_T_US_039_V_4);

    /* Summary */
    printf("\n=== Summary ===\n");
    printf("Tests run: %d, passed: %d, failed: %d\n",
           tests_run, tests_passed, tests_failed);

    return (tests_failed > 0) ? 1 : 0;
}
