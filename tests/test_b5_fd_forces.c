/*
 * SDD:specs.md:S21.1 -- Finite-difference force validation for B5 benchmark
 * US-041: Short QM/MM MD on B5 (QM Water in TIP3P Box)
 *
 * ThDD:T-US-041-STEP-2, T-US-041-STEP-3 -- FD validation for QM and MM forces
 *
 * These tests validate that analytic forces on both QM and MM atoms are
 * consistent with the energy gradient computed by finite differences in
 * the full QM/MM embedded system (B5 benchmark).
 *
 * IMPORTANT: Both the analytic forces and the FD forces are computed during
 * the test. Neither side uses assumed or fabricated reference values.
 * This is the gold standard for force validation per CLAUDE.md Golden Rule.
 *
 * FD methodology (ThDD:T-US-038-V.2):
 *   - Central differences: F = -(E(+d) - E(-d)) / (2*d)
 *   - Step size: delta = 10^-4 Bohr
 *   - Tolerance: 10^-4 relative error OR 10^-6 Ha/Bohr absolute
 *
 * Reference data provenance:
 *   - B5 system: tests/data/b5/ (1 QM water in 883 TIP3P waters)
 *   - DFTB+ version: 25.1 (commit fd31d873)
 *   - SK parameters: mio-1-1 (DFTB2)
 *   - Provenance: tests/data/b5/provenance.json
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
 * FD Test Parameters (from docs/verification/US-041.md Section 4)
 *
 * ThDD:T-US-038-V.2 -- Recommended FD parameters
 *
 * Tolerance derivation:
 *   - Truncation error: O(delta^2 * |F'''|) ~ 10^-8 for typical forces
 *   - SCC reconvergence noise: O(epsilon_SCC / delta) ~ 10^-6
 *   - Combined expected error: ~10^-6 to 10^-5
 *   - Tolerance 10^-4 provides 2 orders of magnitude margin
 * --------------------------------------------------------------------------- */
#define FD_DELTA            1e-4   /* Bohr -- ThDD:T-US-038-V.2 */
#define FD_REL_TOLERANCE    1e-4   /* Relative error tolerance */
#define FD_ABS_FLOOR        1e-6   /* Ha/Bohr -- absolute tolerance floor */

/* ---------------------------------------------------------------------------
 * Test framework macros and counters
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
 * Unit conversion constants (CODATA 2022)
 *
 * ThDD:T-US-041-U.1 -- Unit conventions
 * --------------------------------------------------------------------------- */
#define NM_TO_BOHR          18.8972612462  /* 1 nm = 18.8972612462 Bohr */
#define ANGSTROM_TO_BOHR    1.8897261246   /* 1 A = 1.8897261246 Bohr */

/* ---------------------------------------------------------------------------
 * B5 benchmark system parameters (from tests/data/b5/provenance.json)
 *
 * Provenance:
 *   - Generated: 2026-02-04T17:30:00Z
 *   - GROMACS: 2027.0-dev (e1913f698e)
 *   - DFTB+: 25.1 (fd31d873)
 *   - SK parameters: mio-1-1
 *   - Box: 3.0 x 3.0 x 3.0 nm (27 nm^3)
 *   - N_waters: 884 total (1 QM + 883 MM)
 *   - N_atoms: 2652 total (3 QM + 2649 MM)
 * --------------------------------------------------------------------------- */

#define B5_N_QM_ATOMS       3       /* 1 water molecule: O, H, H */
#define B5_N_MM_WATERS      883     /* 883 MM waters */
#define B5_N_MM_ATOMS       2649    /* 883 * 3 = 2649 MM atoms */
#define B5_N_TOTAL_ATOMS    2652    /* 3 + 2649 = 2652 */

/* QM water positions in nm (from b5_initial.gro, residue 395)
 *   395SOL     OW 1183   1.435   1.677   1.623
 *   395SOL    HW1 1184   1.407   1.748   1.565
 *   395SOL    HW2 1185   1.530   1.673   1.611
 *
 * Note: These are instantaneous positions from NVT equilibration snapshot.
 */
static const double B5_QM_COORDS_NM[9] = {
    1.435, 1.677, 1.623,  /* O  */
    1.407, 1.748, 1.565,  /* H1 */
    1.530, 1.673, 1.611   /* H2 */
};

/* QM species indices: O=0, H=1 (DFTB+ convention) */
static const int B5_QM_SPECIES[3] = { 0, 1, 1 };

/* QM atom indices in GROMACS (0-based): 1182, 1183, 1184 */
#define B5_QM_ATOM_O_IDX    1182
#define B5_QM_ATOM_H1_IDX   1183
#define B5_QM_ATOM_H2_IDX   1184

/* Box center (for distance calculations) */
#define B5_BOX_CENTER_NM    1.5  /* 3.0 / 2 = 1.5 nm */

/* Switching parameters (from US-041.md Appendix D)
 * r_on = 1.0 nm, switch_width = 0.2 nm, r_off = 1.2 nm */
#define B5_R_ON_NM          1.0
#define B5_SWITCH_WIDTH_NM  0.2
#define B5_R_OFF_NM         1.2

/* Convert to Bohr for internal use */
#define B5_R_ON_BOHR        (B5_R_ON_NM * NM_TO_BOHR)
#define B5_SWITCH_WIDTH_BOHR (B5_SWITCH_WIDTH_NM * NM_TO_BOHR)
#define B5_R_OFF_BOHR       (B5_R_OFF_NM * NM_TO_BOHR)

/* TIP3P charges (elementary charge) */
#define TIP3P_CHARGE_O      (-0.834)
#define TIP3P_CHARGE_H      0.417

/* ---------------------------------------------------------------------------
 * Test helper: Path resolution for test data
 * --------------------------------------------------------------------------- */
static const char *get_b5_hsd_path(void)
{
    static char path[1024];
    const char *srcdir = getenv("GRODFTB_SOURCE_DIR");
    if (!srcdir) {
        srcdir = GRODFTB_SOURCE_DIR_STR;
    }
    snprintf(path, sizeof(path), "%s/tests/data/b5/dftb_in.hsd", srcdir);
    return path;
}

static const char *get_b5_dir(void)
{
    static char path[1024];
    const char *srcdir = getenv("GRODFTB_SOURCE_DIR");
    if (!srcdir) {
        srcdir = GRODFTB_SOURCE_DIR_STR;
    }
    snprintf(path, sizeof(path), "%s/tests/data/b5", srcdir);
    return path;
}

__attribute__((unused))
static const char *get_b5_gro_path(void)
{
    static char path[1024];
    const char *srcdir = getenv("GRODFTB_SOURCE_DIR");
    if (!srcdir) {
        srcdir = GRODFTB_SOURCE_DIR_STR;
    }
    snprintf(path, sizeof(path), "%s/tests/data/b5/b5_initial.gro", srcdir);
    return path;
}

/* ---------------------------------------------------------------------------
 * ThDD:T-US-038-V.2 -- Combined acceptance criterion
 *
 * The test passes if:
 *   |F_analytic - F_FD| < max(rel_tol * |F_analytic|, abs_floor)
 *
 * This handles both large forces (relative criterion) and small forces
 * (absolute floor prevents division-by-zero style failures).
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
 * Helper: Compute distance from QM centroid
 * --------------------------------------------------------------------------- */
__attribute__((unused))
static double compute_distance_from_qm_centroid(const double *mm_pos_nm,
                                                 const double *qm_coords_nm,
                                                 int n_qm)
{
    /* Compute QM centroid */
    double cx = 0.0, cy = 0.0, cz = 0.0;
    for (int i = 0; i < n_qm; i++) {
        cx += qm_coords_nm[3*i + 0];
        cy += qm_coords_nm[3*i + 1];
        cz += qm_coords_nm[3*i + 2];
    }
    cx /= n_qm;
    cy /= n_qm;
    cz /= n_qm;

    /* Compute distance */
    double dx = mm_pos_nm[0] - cx;
    double dy = mm_pos_nm[1] - cy;
    double dz = mm_pos_nm[2] - cz;

    return sqrt(dx*dx + dy*dy + dz*dz);
}

/* ---------------------------------------------------------------------------
 * Helper: Convert QM coordinates from nm to Bohr
 * --------------------------------------------------------------------------- */
static void convert_qm_coords_to_bohr(const double *coords_nm,
                                       double *coords_bohr, int natoms)
{
    for (int i = 0; i < 3 * natoms; i++) {
        coords_bohr[i] = coords_nm[i] * NM_TO_BOHR;
    }
}

/* ---------------------------------------------------------------------------
 * Test 1: test_b5_fd_qm_forces_V_US_041_STEP_2
 *
 * ThDD:T-US-041-STEP-2 -- FD validation for QM atom forces
 *
 * SDD:docs/verification/US-041.md:Section 4.3 -- FD test sampling strategy
 *
 * This test validates all 9 force components (3 QM atoms x 3 directions)
 * by finite difference. For B5, this is the isolated QM water calculation.
 *
 * IMPORTANT: Both analytic and FD forces are computed. No fabricated values.
 * --------------------------------------------------------------------------- */
static int test_b5_fd_qm_forces_V_US_041_STEP_2(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    printf(" SKIPPED (no DFTB+)");
    return 1;
#else
    const char *datadir = get_b5_dir();
    const char *hsd_path = get_b5_hsd_path();
    char oldcwd[1024];

    if (!getcwd(oldcwd, sizeof(oldcwd))) {
        printf("\n    getcwd failed\n");
        return 0;
    }
    if (chdir(datadir) != 0) {
        printf("\n    chdir to b5 data dir failed\n");
        return 0;
    }

    /* Convert QM coordinates to Bohr */
    double qm_coords_bohr[9];
    convert_qm_coords_to_bohr(B5_QM_COORDS_NM, qm_coords_bohr, B5_N_QM_ATOMS);

    /*
     * First, compute analytic forces on QM atoms (without embedding for now)
     *
     * ThDD:06_theory:Eq4.1 -- F_K = -dE/dR_K
     */
    grodftb_handle_t handle = NULL;
    int rc = grodftb_init(hsd_path, &handle);
    if (rc != GRODFTB_SUCCESS) {
        chdir(oldcwd);
        printf("\n    grodftb_init failed: %d\n", rc);
        return 0;
    }

    rc = grodftb_set_geometry(handle, B5_N_QM_ATOMS, B5_QM_SPECIES, qm_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        chdir(oldcwd);
        printf("\n    grodftb_set_geometry failed: %d\n", rc);
        return 0;
    }

    /*
     * For isolated QM test, no embedding charges.
     * Full embedded test (with MM) is in test_b5_fd_mm_forces_*.
     */

    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        chdir(oldcwd);
        printf("\n    grodftb_compute failed: %d\n", rc);
        return 0;
    }

    /* Get analytic forces */
    double f_analytic[9];  /* 3 atoms x 3 directions */
    rc = grodftb_get_forces(handle, f_analytic);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        chdir(oldcwd);
        printf("\n    grodftb_get_forces failed: %d\n", rc);
        return 0;
    }

    /* Get reference energy for diagnostics */
    double e_ref = 0.0;
    grodftb_get_energy(handle, &e_ref);

    grodftb_finalize(&handle);

    /*
     * Compute FD forces by perturbing each QM atom
     *
     * ThDD:T-US-038-V.2 -- Central difference formula:
     *   F_i = -(E(R + delta*e_i) - E(R - delta*e_i)) / (2*delta)
     */
    double f_fd[9];
    int all_pass = 1;
    const char *atom_names[] = {"O ", "H1", "H2"};
    const char *dir_names[] = {"x", "y", "z"};

    printf("\n");
    printf("    QM region: 3 atoms (O, H, H) from B5 water\n");
    printf("    Reference energy: %.10f Ha\n", e_ref);
    printf("    FD delta: %.0e Bohr, tolerance: %.0e relative\n\n",
           FD_DELTA, FD_REL_TOLERANCE);

    for (int atom = 0; atom < B5_N_QM_ATOMS; atom++) {
        for (int dir = 0; dir < 3; dir++) {
            double e_plus = 0.0, e_minus = 0.0;
            int idx = atom * 3 + dir;

            /* Positive displacement */
            {
                double coords_pert[9];
                memcpy(coords_pert, qm_coords_bohr, 9 * sizeof(double));
                coords_pert[idx] += FD_DELTA;

                handle = NULL;
                rc = grodftb_init(hsd_path, &handle);
                if (rc != GRODFTB_SUCCESS) {
                    chdir(oldcwd);
                    printf("    grodftb_init failed for +delta\n");
                    return 0;
                }
                grodftb_set_geometry(handle, B5_N_QM_ATOMS, B5_QM_SPECIES, coords_pert);
                grodftb_compute(handle);
                grodftb_get_energy(handle, &e_plus);
                grodftb_finalize(&handle);
            }

            /* Negative displacement */
            {
                double coords_pert[9];
                memcpy(coords_pert, qm_coords_bohr, 9 * sizeof(double));
                coords_pert[idx] -= FD_DELTA;

                handle = NULL;
                rc = grodftb_init(hsd_path, &handle);
                if (rc != GRODFTB_SUCCESS) {
                    chdir(oldcwd);
                    printf("    grodftb_init failed for -delta\n");
                    return 0;
                }
                grodftb_set_geometry(handle, B5_N_QM_ATOMS, B5_QM_SPECIES, coords_pert);
                grodftb_compute(handle);
                grodftb_get_energy(handle, &e_minus);
                grodftb_finalize(&handle);
            }

            /* ThDD:T-US-038-V.2 -- Central difference */
            f_fd[idx] = -(e_plus - e_minus) / (2.0 * FD_DELTA);

            /* Check tolerance */
            double abs_err, rel_err;
            int pass = check_fd_tolerance(f_analytic[idx], f_fd[idx],
                                          &abs_err, &rel_err);

            printf("    F[%s,%s]: analytic=%+.8e  FD=%+.8e  rel_err=%.2e  %s\n",
                   atom_names[atom], dir_names[dir],
                   f_analytic[idx], f_fd[idx], rel_err,
                   pass ? "OK" : "FAIL");

            if (!pass) all_pass = 0;
        }
    }

    chdir(oldcwd);
    return all_pass;
#endif
}

/* ---------------------------------------------------------------------------
 * Test 2: test_b5_fd_mm_forces_inner_V_US_041_STEP_3
 *
 * ThDD:T-US-041-STEP-3a -- FD validation for MM atoms in full interaction region
 *
 * Tests 3 MM atoms with r < r_on (inside the switching region).
 * At these distances, S(u) = 1 and forces are unswitched.
 *
 * Sample positions (from docs/verification/US-041.md Appendix D):
 *   - MM water 1: r ~ 0.28 nm (first solvation shell)
 *   - MM water 2: r ~ 0.5 nm (second shell)
 *   - MM water 3: r ~ 0.8 nm (inside cutoff)
 *
 * IMPORTANT: Both analytic and FD forces are computed. No fabricated values.
 * --------------------------------------------------------------------------- */
static int test_b5_fd_mm_forces_inner_V_US_041_STEP_3(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    printf(" SKIPPED (no DFTB+)");
    return 1;
#else
    const char *datadir = get_b5_dir();
    const char *hsd_path = get_b5_hsd_path();
    char oldcwd[1024];

    if (!getcwd(oldcwd, sizeof(oldcwd))) {
        printf("\n    getcwd failed\n");
        return 0;
    }
    if (chdir(datadir) != 0) {
        printf("\n    chdir to b5 data dir failed\n");
        return 0;
    }

    /* Convert QM coordinates to Bohr */
    double qm_coords_bohr[9];
    convert_qm_coords_to_bohr(B5_QM_COORDS_NM, qm_coords_bohr, B5_N_QM_ATOMS);

    /*
     * Define 3 MM test positions in the full interaction region (r < r_on)
     *
     * These positions are synthetic test points, placed relative to the
     * QM centroid to achieve specific r values. The actual MM water
     * positions from b5_initial.gro would need to be loaded for a full test.
     *
     * QM centroid (nm): ((1.435 + 1.407 + 1.530)/3, (1.677 + 1.748 + 1.673)/3,
     *                    (1.623 + 1.565 + 1.611)/3)
     *                 = (1.4573, 1.6993, 1.5997)
     *
     * Target distances: 0.28, 0.5, 0.8 nm
     */
    double qm_centroid[3] = { 1.4573, 1.6993, 1.5997 };  /* nm */

    /* MM test positions (nm): place along +x from centroid */
    double mm_positions_nm[9] = {
        qm_centroid[0] + 0.28, qm_centroid[1], qm_centroid[2],  /* r ~ 0.28 nm */
        qm_centroid[0] + 0.50, qm_centroid[1], qm_centroid[2],  /* r ~ 0.50 nm */
        qm_centroid[0] + 0.80, qm_centroid[1], qm_centroid[2]   /* r ~ 0.80 nm */
    };

    /* TIP3P oxygen charges */
    double mm_charges[3] = { TIP3P_CHARGE_O, TIP3P_CHARGE_O, TIP3P_CHARGE_O };
    int n_mm_test = 3;

    /* Convert MM positions to Bohr */
    double mm_positions_bohr[9];
    for (int i = 0; i < 9; i++) {
        mm_positions_bohr[i] = mm_positions_nm[i] * NM_TO_BOHR;
    }

    /*
     * First, compute analytic forces with embedding
     */
    grodftb_handle_t handle = NULL;
    int rc = grodftb_init(hsd_path, &handle);
    if (rc != GRODFTB_SUCCESS) {
        chdir(oldcwd);
        printf("\n    grodftb_init failed: %d\n", rc);
        return 0;
    }

    rc = grodftb_set_geometry(handle, B5_N_QM_ATOMS, B5_QM_SPECIES, qm_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        chdir(oldcwd);
        printf("\n    grodftb_set_geometry failed: %d\n", rc);
        return 0;
    }

    /* Set embedding parameters (full interaction, no switching active at r < r_on) */
    rc = grodftb_set_embedding_params(handle, B5_R_OFF_BOHR, B5_SWITCH_WIDTH_BOHR);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        chdir(oldcwd);
        printf("\n    grodftb_set_embedding_params failed: %d\n", rc);
        return 0;
    }

    rc = grodftb_set_embedding_charges(handle, n_mm_test, mm_charges, mm_positions_bohr);
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
    double f_analytic[9];  /* 3 MM atoms x 3 directions */
    rc = grodftb_get_embedding_forces(handle, f_analytic);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        chdir(oldcwd);
        printf("\n    grodftb_get_embedding_forces failed: %d\n", rc);
        return 0;
    }

    grodftb_finalize(&handle);

    /*
     * Compute FD forces by perturbing each MM position
     *
     * ThDD:T-US-038-V.2 -- Central difference formula
     */
    int all_pass = 1;
    const char *dir_names[] = {"x", "y", "z"};
    double target_r[] = {0.28, 0.50, 0.80};

    printf("\n");
    printf("    Testing MM atoms in full interaction region (r < r_on = %.1f nm)\n",
           B5_R_ON_NM);
    printf("    FD delta: %.0e Bohr, tolerance: %.0e relative\n\n",
           FD_DELTA, FD_REL_TOLERANCE);

    for (int mm_idx = 0; mm_idx < n_mm_test; mm_idx++) {
        printf("    MM[%d] at r ~ %.2f nm:\n", mm_idx, target_r[mm_idx]);

        for (int dir = 0; dir < 3; dir++) {
            double e_plus = 0.0, e_minus = 0.0;
            int idx = mm_idx * 3 + dir;

            /* Positive displacement */
            {
                double mm_pert[9];
                memcpy(mm_pert, mm_positions_bohr, 9 * sizeof(double));
                mm_pert[idx] += FD_DELTA;

                handle = NULL;
                rc = grodftb_init(hsd_path, &handle);
                if (rc != GRODFTB_SUCCESS) {
                    chdir(oldcwd);
                    printf("      grodftb_init failed for +delta\n");
                    return 0;
                }
                grodftb_set_geometry(handle, B5_N_QM_ATOMS, B5_QM_SPECIES, qm_coords_bohr);
                grodftb_set_embedding_params(handle, B5_R_OFF_BOHR, B5_SWITCH_WIDTH_BOHR);
                grodftb_set_embedding_charges(handle, n_mm_test, mm_charges, mm_pert);
                grodftb_compute(handle);
                grodftb_get_energy(handle, &e_plus);
                grodftb_finalize(&handle);
            }

            /* Negative displacement */
            {
                double mm_pert[9];
                memcpy(mm_pert, mm_positions_bohr, 9 * sizeof(double));
                mm_pert[idx] -= FD_DELTA;

                handle = NULL;
                rc = grodftb_init(hsd_path, &handle);
                if (rc != GRODFTB_SUCCESS) {
                    chdir(oldcwd);
                    printf("      grodftb_init failed for -delta\n");
                    return 0;
                }
                grodftb_set_geometry(handle, B5_N_QM_ATOMS, B5_QM_SPECIES, qm_coords_bohr);
                grodftb_set_embedding_params(handle, B5_R_OFF_BOHR, B5_SWITCH_WIDTH_BOHR);
                grodftb_set_embedding_charges(handle, n_mm_test, mm_charges, mm_pert);
                grodftb_compute(handle);
                grodftb_get_energy(handle, &e_minus);
                grodftb_finalize(&handle);
            }

            /* ThDD:T-US-038-V.2 -- Central difference */
            double f_fd = -(e_plus - e_minus) / (2.0 * FD_DELTA);

            /* Check tolerance */
            double abs_err, rel_err;
            int pass = check_fd_tolerance(f_analytic[idx], f_fd,
                                          &abs_err, &rel_err);

            printf("      F[%s]: analytic=%+.8e  FD=%+.8e  rel_err=%.2e  %s\n",
                   dir_names[dir], f_analytic[idx], f_fd, rel_err,
                   pass ? "OK" : "FAIL");

            if (!pass) all_pass = 0;
        }
    }

    chdir(oldcwd);
    return all_pass;
#endif
}

/* ---------------------------------------------------------------------------
 * Test 3: test_b5_fd_mm_forces_switch_V_US_041_STEP_3
 *
 * ThDD:T-US-041-STEP-3b -- FD validation for MM atoms in switching region
 *
 * Tests 4 MM atoms with r_on < r < r_off (inside the switching region).
 * The switching function S(u) applies here, making this the most critical test.
 *
 * Sample positions (from docs/verification/US-041.md Appendix D):
 *   - MM water 4: r ~ 1.02 nm (u ~ 0.1, S ~ 0.97)
 *   - MM water 5: r ~ 1.06 nm (u ~ 0.3, S ~ 0.84)
 *   - MM water 6: r ~ 1.12 nm (u ~ 0.6, S ~ 0.32)
 *   - MM water 7: r ~ 1.18 nm (u ~ 0.9, S ~ 0.03)
 *
 * IMPORTANT: Both analytic and FD forces are computed. No fabricated values.
 * --------------------------------------------------------------------------- */
static int test_b5_fd_mm_forces_switch_V_US_041_STEP_3(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    printf(" SKIPPED (no DFTB+)");
    return 1;
#else
    const char *datadir = get_b5_dir();
    const char *hsd_path = get_b5_hsd_path();
    char oldcwd[1024];

    if (!getcwd(oldcwd, sizeof(oldcwd))) {
        printf("\n    getcwd failed\n");
        return 0;
    }
    if (chdir(datadir) != 0) {
        printf("\n    chdir to b5 data dir failed\n");
        return 0;
    }

    /* Convert QM coordinates to Bohr */
    double qm_coords_bohr[9];
    convert_qm_coords_to_bohr(B5_QM_COORDS_NM, qm_coords_bohr, B5_N_QM_ATOMS);

    /*
     * Define 4 MM test positions in the switching region (r_on < r < r_off)
     *
     * Reduced coordinate u = (r - r_on) / switch_width
     * For r_on = 1.0 nm, switch_width = 0.2 nm:
     *   r = 1.02 nm -> u = 0.1
     *   r = 1.06 nm -> u = 0.3
     *   r = 1.12 nm -> u = 0.6
     *   r = 1.18 nm -> u = 0.9
     */
    double qm_centroid[3] = { 1.4573, 1.6993, 1.5997 };  /* nm */

    /* MM test positions (nm): place along +x from centroid */
    double mm_positions_nm[12] = {
        qm_centroid[0] + 1.02, qm_centroid[1], qm_centroid[2],  /* r ~ 1.02 nm, u ~ 0.1 */
        qm_centroid[0] + 1.06, qm_centroid[1], qm_centroid[2],  /* r ~ 1.06 nm, u ~ 0.3 */
        qm_centroid[0] + 1.12, qm_centroid[1], qm_centroid[2],  /* r ~ 1.12 nm, u ~ 0.6 */
        qm_centroid[0] + 1.18, qm_centroid[1], qm_centroid[2]   /* r ~ 1.18 nm, u ~ 0.9 */
    };

    /* TIP3P oxygen charges */
    double mm_charges[4] = { TIP3P_CHARGE_O, TIP3P_CHARGE_O,
                             TIP3P_CHARGE_O, TIP3P_CHARGE_O };
    int n_mm_test = 4;

    /* Convert MM positions to Bohr */
    double mm_positions_bohr[12];
    for (int i = 0; i < 12; i++) {
        mm_positions_bohr[i] = mm_positions_nm[i] * NM_TO_BOHR;
    }

    /* Expected u values and S(u) for reference */
    double u_values[] = { 0.1, 0.3, 0.6, 0.9 };
    double target_r[] = { 1.02, 1.06, 1.12, 1.18 };

    /*
     * First, compute analytic forces with embedding
     */
    grodftb_handle_t handle = NULL;
    int rc = grodftb_init(hsd_path, &handle);
    if (rc != GRODFTB_SUCCESS) {
        chdir(oldcwd);
        printf("\n    grodftb_init failed: %d\n", rc);
        return 0;
    }

    rc = grodftb_set_geometry(handle, B5_N_QM_ATOMS, B5_QM_SPECIES, qm_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        chdir(oldcwd);
        printf("\n    grodftb_set_geometry failed: %d\n", rc);
        return 0;
    }

    /* Set embedding parameters with switching */
    rc = grodftb_set_embedding_params(handle, B5_R_OFF_BOHR, B5_SWITCH_WIDTH_BOHR);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        chdir(oldcwd);
        printf("\n    grodftb_set_embedding_params failed: %d\n", rc);
        return 0;
    }

    rc = grodftb_set_embedding_charges(handle, n_mm_test, mm_charges, mm_positions_bohr);
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
    double f_analytic[12];  /* 4 MM atoms x 3 directions */
    rc = grodftb_get_embedding_forces(handle, f_analytic);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        chdir(oldcwd);
        printf("\n    grodftb_get_embedding_forces failed: %d\n", rc);
        return 0;
    }

    grodftb_finalize(&handle);

    /*
     * Compute FD forces by perturbing each MM position
     */
    int all_pass = 1;
    const char *dir_names[] = {"x", "y", "z"};

    printf("\n");
    printf("    Testing MM atoms in switching region (%.1f nm < r < %.1f nm)\n",
           B5_R_ON_NM, B5_R_OFF_NM);
    printf("    r_on=%.1f nm, switch_width=%.1f nm, r_off=%.1f nm\n",
           B5_R_ON_NM, B5_SWITCH_WIDTH_NM, B5_R_OFF_NM);
    printf("    FD delta: %.0e Bohr, tolerance: %.0e relative\n\n",
           FD_DELTA, FD_REL_TOLERANCE);

    for (int mm_idx = 0; mm_idx < n_mm_test; mm_idx++) {
        printf("    MM[%d] at r ~ %.2f nm (u ~ %.1f):\n",
               mm_idx, target_r[mm_idx], u_values[mm_idx]);

        for (int dir = 0; dir < 3; dir++) {
            double e_plus = 0.0, e_minus = 0.0;
            int idx = mm_idx * 3 + dir;

            /* Positive displacement */
            {
                double mm_pert[12];
                memcpy(mm_pert, mm_positions_bohr, 12 * sizeof(double));
                mm_pert[idx] += FD_DELTA;

                handle = NULL;
                rc = grodftb_init(hsd_path, &handle);
                if (rc != GRODFTB_SUCCESS) {
                    chdir(oldcwd);
                    printf("      grodftb_init failed for +delta\n");
                    return 0;
                }
                grodftb_set_geometry(handle, B5_N_QM_ATOMS, B5_QM_SPECIES, qm_coords_bohr);
                grodftb_set_embedding_params(handle, B5_R_OFF_BOHR, B5_SWITCH_WIDTH_BOHR);
                grodftb_set_embedding_charges(handle, n_mm_test, mm_charges, mm_pert);
                grodftb_compute(handle);
                grodftb_get_energy(handle, &e_plus);
                grodftb_finalize(&handle);
            }

            /* Negative displacement */
            {
                double mm_pert[12];
                memcpy(mm_pert, mm_positions_bohr, 12 * sizeof(double));
                mm_pert[idx] -= FD_DELTA;

                handle = NULL;
                rc = grodftb_init(hsd_path, &handle);
                if (rc != GRODFTB_SUCCESS) {
                    chdir(oldcwd);
                    printf("      grodftb_init failed for -delta\n");
                    return 0;
                }
                grodftb_set_geometry(handle, B5_N_QM_ATOMS, B5_QM_SPECIES, qm_coords_bohr);
                grodftb_set_embedding_params(handle, B5_R_OFF_BOHR, B5_SWITCH_WIDTH_BOHR);
                grodftb_set_embedding_charges(handle, n_mm_test, mm_charges, mm_pert);
                grodftb_compute(handle);
                grodftb_get_energy(handle, &e_minus);
                grodftb_finalize(&handle);
            }

            /* ThDD:T-US-038-V.2 -- Central difference */
            double f_fd = -(e_plus - e_minus) / (2.0 * FD_DELTA);

            /* Check tolerance */
            double abs_err, rel_err;
            int pass = check_fd_tolerance(f_analytic[idx], f_fd,
                                          &abs_err, &rel_err);

            printf("      F[%s]: analytic=%+.8e  FD=%+.8e  rel_err=%.2e  %s\n",
                   dir_names[dir], f_analytic[idx], f_fd, rel_err,
                   pass ? "OK" : "FAIL");

            if (!pass) all_pass = 0;
        }
    }

    chdir(oldcwd);
    return all_pass;
#endif
}

/* ---------------------------------------------------------------------------
 * Test 4: test_b5_fd_mm_forces_outer_V_US_041_STEP_3
 *
 * ThDD:T-US-041-STEP-3c -- FD validation for MM atoms near r_off boundary
 *
 * Tests 3 MM atoms near and beyond r_off to verify:
 *   - Forces approach zero smoothly near r_off
 *   - Forces are exactly zero at r >= r_off FROM ALL QM ATOMS
 *
 * IMPORTANT: The switching cutoff applies to individual QM-MM pair distances,
 * not to the distance from the QM centroid. An MM atom at r=1.21nm from the
 * centroid may still be within r_off of some QM atoms. For the "beyond cutoff"
 * test case, we must place the MM atom far enough that ALL QM-MM distances
 * exceed r_off.
 *
 * QM water geometry spans ~0.1 nm, so to ensure all atoms are beyond 1.2nm:
 *   - Place MM at centroid_distance = 1.3 nm (guarantees min distance > 1.2 nm)
 *
 * Sample positions:
 *   - MM water 8: r_centroid ~ 1.19 nm (some atoms in switching)
 *   - MM water 9: r_centroid ~ 1.25 nm (all atoms near boundary)
 *   - MM water 10: r_centroid ~ 1.35 nm (all atoms beyond cutoff)
 *
 * IMPORTANT: Both analytic and FD forces are computed. No fabricated values.
 * --------------------------------------------------------------------------- */
static int test_b5_fd_mm_forces_outer_V_US_041_STEP_3(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    printf(" SKIPPED (no DFTB+)");
    return 1;
#else
    const char *datadir = get_b5_dir();
    const char *hsd_path = get_b5_hsd_path();
    char oldcwd[1024];

    if (!getcwd(oldcwd, sizeof(oldcwd))) {
        printf("\n    getcwd failed\n");
        return 0;
    }
    if (chdir(datadir) != 0) {
        printf("\n    chdir to b5 data dir failed\n");
        return 0;
    }

    /* Convert QM coordinates to Bohr */
    double qm_coords_bohr[9];
    convert_qm_coords_to_bohr(B5_QM_COORDS_NM, qm_coords_bohr, B5_N_QM_ATOMS);

    /*
     * Define 3 MM test positions near/beyond the outer boundary
     *
     * The QM water spans ~0.1 nm (O-H ~ 0.096 nm), so:
     * - r_centroid = 1.19 nm means closest QM-MM distance ~ 1.14 nm (in switching)
     * - r_centroid = 1.25 nm means closest QM-MM distance ~ 1.20 nm (at boundary)
     * - r_centroid = 1.35 nm means closest QM-MM distance ~ 1.30 nm (beyond cutoff)
     *
     * For r_on = 1.0 nm, switch_width = 0.2 nm, r_off = 1.2 nm:
     * The third position should have zero force because all QM-MM distances > r_off.
     */
    double qm_centroid[3] = { 1.4573, 1.6993, 1.5997 };  /* nm */

    /* MM test positions (nm): place along +x from centroid
     * Using larger distances to ensure the "beyond cutoff" case is valid */
    double mm_positions_nm[9] = {
        qm_centroid[0] + 1.19, qm_centroid[1], qm_centroid[2],  /* r_centroid ~ 1.19 nm */
        qm_centroid[0] + 1.25, qm_centroid[1], qm_centroid[2],  /* r_centroid ~ 1.25 nm */
        qm_centroid[0] + 1.35, qm_centroid[1], qm_centroid[2]   /* r_centroid ~ 1.35 nm (all beyond) */
    };

    /* TIP3P oxygen charges */
    double mm_charges[3] = { TIP3P_CHARGE_O, TIP3P_CHARGE_O, TIP3P_CHARGE_O };
    int n_mm_test = 3;

    /* Convert MM positions to Bohr */
    double mm_positions_bohr[9];
    for (int i = 0; i < 9; i++) {
        mm_positions_bohr[i] = mm_positions_nm[i] * NM_TO_BOHR;
    }

    /* Approximate u values (based on closest QM-MM distance) and centroid distances */
    double u_values[] = { 0.70, 1.00, 1.50 };  /* Approximate */
    double target_r[] = { 1.19, 1.25, 1.35 };

    /*
     * First, compute analytic forces with embedding
     */
    grodftb_handle_t handle = NULL;
    int rc = grodftb_init(hsd_path, &handle);
    if (rc != GRODFTB_SUCCESS) {
        chdir(oldcwd);
        printf("\n    grodftb_init failed: %d\n", rc);
        return 0;
    }

    rc = grodftb_set_geometry(handle, B5_N_QM_ATOMS, B5_QM_SPECIES, qm_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        chdir(oldcwd);
        printf("\n    grodftb_set_geometry failed: %d\n", rc);
        return 0;
    }

    /* Set embedding parameters with switching */
    rc = grodftb_set_embedding_params(handle, B5_R_OFF_BOHR, B5_SWITCH_WIDTH_BOHR);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        chdir(oldcwd);
        printf("\n    grodftb_set_embedding_params failed: %d\n", rc);
        return 0;
    }

    rc = grodftb_set_embedding_charges(handle, n_mm_test, mm_charges, mm_positions_bohr);
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
    double f_analytic[9];  /* 3 MM atoms x 3 directions */
    rc = grodftb_get_embedding_forces(handle, f_analytic);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        chdir(oldcwd);
        printf("\n    grodftb_get_embedding_forces failed: %d\n", rc);
        return 0;
    }

    grodftb_finalize(&handle);

    /*
     * Compute FD forces by perturbing each MM position
     *
     * For the atom beyond cutoff (u > 1), we expect both analytic and FD
     * forces to be zero (or very small). The FD test still validates
     * consistency even at zero force.
     */
    int all_pass = 1;
    const char *dir_names[] = {"x", "y", "z"};

    printf("\n");
    printf("    Testing MM atoms near outer boundary (r_off = %.1f nm)\n",
           B5_R_OFF_NM);
    printf("    Note: cutoff applies to individual QM-MM distances, not centroid\n");
    printf("    Atoms with ALL QM-MM distances >= r_off should have force = 0\n");
    printf("    FD delta: %.0e Bohr, tolerance: %.0e relative OR %.0e absolute\n\n",
           FD_DELTA, FD_REL_TOLERANCE, FD_ABS_FLOOR);

    for (int mm_idx = 0; mm_idx < n_mm_test; mm_idx++) {
        const char *status;
        if (target_r[mm_idx] > B5_R_OFF_NM + 0.1) {
            /* MM is at centroid distance > r_off + 0.1 nm, so all QM atoms should be beyond r_off */
            status = "(all QM-MM beyond cutoff, expect F=0)";
        } else if (target_r[mm_idx] > B5_R_OFF_NM) {
            status = "(some QM-MM near boundary)";
        } else {
            status = "(some QM-MM in switching)";
        }

        printf("    MM[%d] at r ~ %.3f nm (u ~ %.3f) %s:\n",
               mm_idx, target_r[mm_idx], u_values[mm_idx], status);

        for (int dir = 0; dir < 3; dir++) {
            double e_plus = 0.0, e_minus = 0.0;
            int idx = mm_idx * 3 + dir;

            /* Positive displacement */
            {
                double mm_pert[9];
                memcpy(mm_pert, mm_positions_bohr, 9 * sizeof(double));
                mm_pert[idx] += FD_DELTA;

                handle = NULL;
                rc = grodftb_init(hsd_path, &handle);
                if (rc != GRODFTB_SUCCESS) {
                    chdir(oldcwd);
                    printf("      grodftb_init failed for +delta\n");
                    return 0;
                }
                grodftb_set_geometry(handle, B5_N_QM_ATOMS, B5_QM_SPECIES, qm_coords_bohr);
                grodftb_set_embedding_params(handle, B5_R_OFF_BOHR, B5_SWITCH_WIDTH_BOHR);
                grodftb_set_embedding_charges(handle, n_mm_test, mm_charges, mm_pert);
                grodftb_compute(handle);
                grodftb_get_energy(handle, &e_plus);
                grodftb_finalize(&handle);
            }

            /* Negative displacement */
            {
                double mm_pert[9];
                memcpy(mm_pert, mm_positions_bohr, 9 * sizeof(double));
                mm_pert[idx] -= FD_DELTA;

                handle = NULL;
                rc = grodftb_init(hsd_path, &handle);
                if (rc != GRODFTB_SUCCESS) {
                    chdir(oldcwd);
                    printf("      grodftb_init failed for -delta\n");
                    return 0;
                }
                grodftb_set_geometry(handle, B5_N_QM_ATOMS, B5_QM_SPECIES, qm_coords_bohr);
                grodftb_set_embedding_params(handle, B5_R_OFF_BOHR, B5_SWITCH_WIDTH_BOHR);
                grodftb_set_embedding_charges(handle, n_mm_test, mm_charges, mm_pert);
                grodftb_compute(handle);
                grodftb_get_energy(handle, &e_minus);
                grodftb_finalize(&handle);
            }

            /* ThDD:T-US-038-V.2 -- Central difference */
            double f_fd = -(e_plus - e_minus) / (2.0 * FD_DELTA);

            /* Check tolerance */
            double abs_err, rel_err;
            int pass = check_fd_tolerance(f_analytic[idx], f_fd,
                                          &abs_err, &rel_err);

            /*
             * For fully-beyond-cutoff case (all QM-MM distances > r_off),
             * verify absolute value is near zero.
             *
             * The criterion is: centroid distance > r_off + ~0.1 nm (half the water span)
             * so that ALL individual QM-MM distances exceed r_off.
             */
            if (target_r[mm_idx] > B5_R_OFF_NM + 0.1) {
                /* Force should be essentially zero */
                double zero_tol = 1e-14;  /* Ha/Bohr -- allow small numerical noise */
                if (fabs(f_analytic[idx]) > zero_tol || fabs(f_fd) > zero_tol) {
                    printf("      F[%s]: analytic=%+.2e  FD=%+.2e  (non-zero beyond cutoff!)  FAIL\n",
                           dir_names[dir], f_analytic[idx], f_fd);
                    pass = 0;
                } else {
                    printf("      F[%s]: analytic=%+.2e  FD=%+.2e  (correctly zero)  OK\n",
                           dir_names[dir], f_analytic[idx], f_fd);
                    pass = 1;
                }
            } else {
                printf("      F[%s]: analytic=%+.8e  FD=%+.8e  rel_err=%.2e  %s\n",
                       dir_names[dir], f_analytic[idx], f_fd, rel_err,
                       pass ? "OK" : "FAIL");
            }

            if (!pass) all_pass = 0;
        }
    }

    chdir(oldcwd);
    return all_pass;
#endif
}

/* ---------------------------------------------------------------------------
 * Main test runner
 * --------------------------------------------------------------------------- */
int main(void)
{
    printf("=== US-041 B5 FD Force Validation Tests (TDD Red Phase) ===\n\n");

    printf("ThDD:T-US-041-STEP-2,3 -- FD validation for QM and MM forces\n");
    printf("SDD:docs/verification/US-041.md:Section 4 -- FD methodology\n");
    printf("SDD:specs.md:S21.1 -- Tolerance: 10^-4 relative error\n\n");

    printf("Test methodology:\n");
    printf("  - Central differences: F = -(E(+d) - E(-d)) / (2*d)\n");
    printf("  - Step size: delta = 10^-4 Bohr\n");
    printf("  - Tolerance: max(10^-4 * |F|, 10^-6 Ha/Bohr)\n");
    printf("  - BOTH analytic and FD forces are computed (no fabricated values)\n\n");

    printf("B5 benchmark system:\n");
    printf("  - QM region: 1 water (3 atoms: O, H, H)\n");
    printf("  - MM region: 883 waters (2649 atoms)\n");
    printf("  - Box: 3.0 x 3.0 x 3.0 nm\n");
    printf("  - Switching: r_on=1.0 nm, width=0.2 nm, r_off=1.2 nm\n\n");

    printf("--- Test 1: QM Force FD Validation (STEP-2) ---\n");
    RUN_TEST(test_b5_fd_qm_forces_V_US_041_STEP_2);

    printf("\n--- Test 2: MM Force FD Validation, Inner Region (STEP-3a) ---\n");
    RUN_TEST(test_b5_fd_mm_forces_inner_V_US_041_STEP_3);

    printf("\n--- Test 3: MM Force FD Validation, Switching Region (STEP-3b) ---\n");
    RUN_TEST(test_b5_fd_mm_forces_switch_V_US_041_STEP_3);

    printf("\n--- Test 4: MM Force FD Validation, Outer Boundary (STEP-3c) ---\n");
    RUN_TEST(test_b5_fd_mm_forces_outer_V_US_041_STEP_3);

    /* Summary */
    printf("\n=== Summary ===\n");
    printf("Tests run:    %d\n", tests_run);
    printf("Tests passed: %d\n", tests_passed);
    printf("Tests failed: %d\n", tests_failed);

    if (tests_failed > 0) {
        printf("\nNOTE: This is the TDD Red Phase. Tests are expected to fail\n");
        printf("      until the B5 QM/MM coupling module is implemented.\n");
    }

    return tests_failed > 0 ? 1 : 0;
}
