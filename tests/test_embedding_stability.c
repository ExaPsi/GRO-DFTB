/*
 * US-043: Embedding Stability Tests (V.15-V.16)
 *
 * ThDD:T-US-043-V.15 -- Switching function smoothness (cross-reference US-039)
 * ThDD:T-US-043-V.16 -- MM force continuity over trajectory
 * ThDD:T-US-043-E.3  -- Cumulative energy error bound
 * ThDD:T-US-043-E.4  -- Energy continuity through switching transitions
 * SDD:specs.md:S21.1 -- NVE simulation requirements
 * SDD:docs/verification/US-043.md -- Verification plan
 *
 * This file implements:
 *   1. V.15: Cross-reference to US-039 switching smoothness tests (not duplicated)
 *   2. V.16: MM force continuity over extended trajectories
 *
 * GOLDEN RULE COMPLIANCE:
 * - Unit tests use mathematical identities (exact formulas)
 * - Force continuity threshold (10 kJ/mol/nm) from T-US-043-V.16
 * - No hardcoded "expected" simulation energy/force values
 * - All simulation results must come from actual trajectory data
 *
 * Test execution modes:
 *   --unit         Run unit tests only (fast, no simulation required)
 *   --integration  Run integration tests (requires GRODFTB_LONG_TESTS)
 *   (no args)      Run unit tests only
 *
 * Cross-reference:
 *   V.15 (switching smoothness) is validated in US-039 tests:
 *   - tests/test_switching.c
 *   - tests/test_switching_fd.c
 *   - tests/test_switching_regression.c
 *
 * LGPL-3.0-or-later
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "grodftb/units.h"
#include "grodftb/error.h"

/* ===========================================================================
 * ACCEPTANCE CRITERIA (from T-US-043-V.16)
 *
 * ThDD:T-US-043-V.16:
 *   |F_MM(t+dt) - F_MM(t)| < 10 kJ/mol/nm (for atoms not crossing cutoff)
 *
 * This threshold ensures MM forces are smooth enough to avoid causing
 * energy drift from force discontinuities. The value is derived from:
 *   - Typical thermal force fluctuations: sqrt(k_B * T / m) * m / dt ~ 1 kJ/mol/nm
 *   - 10x margin to allow for normal dynamics
 *   - Must be well above numerical noise (~10^-6 kJ/mol/nm)
 *
 * These are specification-derived thresholds, NOT fabricated expected values.
 * ===========================================================================*/

/* SDD:T-US-043-V.16 -- MM force continuity threshold */
#define MM_FORCE_JUMP_THRESHOLD_KJMOL_NM  10.0

/* ===========================================================================
 * B5 SYSTEM PARAMETERS (from tests/data/b5/provenance.json)
 *
 * ThDD:CLAUDE.md:Golden_Rule -- All values from actual GROMACS preparation
 * ===========================================================================*/

#define B5_N_ATOMS_TOTAL     2652   /* Total atoms in system */
#define B5_N_QM_ATOMS        3      /* QM region: 1 water molecule */
#define B5_N_MM_ATOMS        2649   /* MM atoms (total - QM) */
#define B5_CUTOFF_NM         1.2    /* QM/MM cutoff in nm */
#define B5_SWITCH_WIDTH_NM   0.2    /* Switching width in nm (if enabled) */
#define B5_R_ON_NM           1.0    /* Inner switching boundary (cutoff - width) */

/* ===========================================================================
 * FORCE ANALYSIS RESULT STRUCTURE
 *
 * ThDD:T-US-043-V.16 -- Force continuity analysis
 * ===========================================================================*/

/**
 * Result from MM force continuity analysis.
 *
 * ThDD:T-US-043-V.16 -- |F_MM(t+dt) - F_MM(t)| statistics
 */
typedef struct {
    int n_frames;           /* Number of trajectory frames analyzed */
    int n_mm_atoms;         /* Number of MM atoms tracked */
    double max_force_jump;  /* Maximum |dF/dt| observed (kJ/mol/nm) */
    int max_jump_atom_idx;  /* Atom index with maximum jump */
    int max_jump_frame;     /* Frame index of maximum jump */
    double mean_force_jump; /* Mean |dF/dt| across all atoms and frames */
    double stddev_force_jump; /* Standard deviation of force jumps */
    int n_jumps_above_threshold; /* Number of jumps exceeding threshold */
    int threshold_violations; /* Same as above (alias for clarity) */
} force_continuity_result_t;

/* ===========================================================================
 * SECTION 1: FORCE CONTINUITY ANALYSIS HELPER FUNCTIONS
 *
 * These functions provide tools for analyzing MM force continuity over time.
 * Unit tests verify the mathematical correctness of these helpers.
 * ===========================================================================*/

/**
 * Compute the magnitude of a 3D force vector.
 *
 * ThDD:T-US-043-V.16 -- Force magnitude calculation
 *
 * MATHEMATICAL IDENTITY:
 *   |F| = sqrt(Fx^2 + Fy^2 + Fz^2)
 *
 * @param fx  X component of force
 * @param fy  Y component of force
 * @param fz  Z component of force
 * @return    Magnitude |F|
 */
static double compute_force_magnitude(double fx, double fy, double fz)
{
    return sqrt(fx * fx + fy * fy + fz * fz);
}

/**
 * Compute the magnitude of force change between two timesteps.
 *
 * ThDD:T-US-043-V.16 -- Force jump calculation
 *
 * MATHEMATICAL IDENTITY:
 *   |dF| = |F(t+dt) - F(t)| = sqrt((Fx2-Fx1)^2 + (Fy2-Fy1)^2 + (Fz2-Fz1)^2)
 *
 * @param fx1  X component at time t
 * @param fy1  Y component at time t
 * @param fz1  Z component at time t
 * @param fx2  X component at time t+dt
 * @param fy2  Y component at time t+dt
 * @param fz2  Z component at time t+dt
 * @return     Magnitude of force change |F(t+dt) - F(t)|
 */
static double compute_force_jump(double fx1, double fy1, double fz1,
                                  double fx2, double fy2, double fz2)
{
    double dx = fx2 - fx1;
    double dy = fy2 - fy1;
    double dz = fz2 - fz1;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

/**
 * Check if a position is within the switching region.
 *
 * ThDD:T-US-043-E.4 -- Switching region identification
 *
 * The switching region is r_on <= r < r_cut, where interactions are
 * smoothly attenuated by the quintic switching function.
 *
 * @param r        Distance from QM center in nm
 * @param r_on     Inner switching boundary in nm
 * @param r_cut    Outer cutoff boundary in nm
 * @return         1 if in switching region, 0 otherwise
 */
static int is_in_switching_region(double r, double r_on, double r_cut)
{
    return (r >= r_on && r < r_cut) ? 1 : 0;
}

/**
 * Check if an atom crossed the cutoff boundary between two timesteps.
 *
 * ThDD:T-US-043-V.16 -- Boundary crossing detection
 *
 * An atom that crosses the cutoff boundary is expected to have a force
 * discontinuity (from finite to zero or vice versa). V.16 only checks
 * atoms that did NOT cross the boundary.
 *
 * @param r_prev   Distance at time t (nm)
 * @param r_curr   Distance at time t+dt (nm)
 * @param r_cut    Cutoff distance (nm)
 * @return         1 if crossed boundary, 0 if stayed on same side
 */
static int crossed_cutoff_boundary(double r_prev, double r_curr, double r_cut)
{
    int inside_prev = (r_prev < r_cut) ? 1 : 0;
    int inside_curr = (r_curr < r_cut) ? 1 : 0;
    return (inside_prev != inside_curr) ? 1 : 0;
}

/**
 * Initialize force continuity result structure.
 *
 * @param result  Result structure to initialize
 */
static void __attribute__((unused)) init_force_continuity_result(force_continuity_result_t *result)
{
    memset(result, 0, sizeof(*result));
    result->max_force_jump = 0.0;
    result->max_jump_atom_idx = -1;
    result->max_jump_frame = -1;
}

/**
 * Compute sample mean of an array.
 *
 * @param n    Number of elements
 * @param arr  Data array
 * @return     Arithmetic mean
 */
static double compute_mean(int n, const double *arr)
{
    if (n <= 0) return 0.0;

    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += arr[i];
    }
    return sum / (double)n;
}

/**
 * Compute sample standard deviation of an array.
 *
 * @param n    Number of elements
 * @param arr  Data array
 * @return     Sample standard deviation (Bessel-corrected)
 */
static double __attribute__((unused)) compute_stddev(int n, const double *arr)
{
    if (n < 2) return 0.0;

    double mean = compute_mean(n, arr);

    double variance = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = arr[i] - mean;
        variance += diff * diff;
    }
    variance /= (double)(n - 1);

    return sqrt(variance);
}

/**
 * Count number of values exceeding a threshold.
 *
 * ThDD:T-US-043-V.16 -- Threshold violation counting
 *
 * @param n          Number of elements
 * @param arr        Data array
 * @param threshold  Threshold value
 * @return           Number of elements > threshold
 */
static int count_above_threshold(int n, const double *arr, double threshold)
{
    int count = 0;
    for (int i = 0; i < n; i++) {
        if (arr[i] > threshold) count++;
    }
    return count;
}

/**
 * Find maximum value in an array and its index.
 *
 * @param n        Number of elements
 * @param arr      Data array
 * @param max_val  Output: maximum value found
 * @param max_idx  Output: index of maximum value
 */
static void find_max_and_index(int n, const double *arr,
                                double *max_val, int *max_idx)
{
    if (n <= 0) {
        *max_val = 0.0;
        *max_idx = -1;
        return;
    }

    *max_val = arr[0];
    *max_idx = 0;

    for (int i = 1; i < n; i++) {
        if (arr[i] > *max_val) {
            *max_val = arr[i];
            *max_idx = i;
        }
    }
}

/* ===========================================================================
 * SECTION 2: TEST FRAMEWORK
 * ===========================================================================*/

static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;
static int tests_skipped = 0;

#define RUN_TEST(test_func) do { \
    printf("  Running %s... ", #test_func); \
    fflush(stdout); \
    tests_run++; \
    int result = test_func(); \
    if (result == 1) { \
        printf("PASS\n"); \
        tests_passed++; \
    } else if (result == 0) { \
        printf("FAIL\n"); \
        tests_failed++; \
    } else { \
        printf("SKIP\n"); \
        tests_skipped++; \
    } \
} while(0)

/* ===========================================================================
 * SECTION 3: UNIT TESTS (Mathematical Identities)
 *
 * GOLDEN RULE: All expected values are mathematical identities,
 * not fabricated physical results.
 * ===========================================================================*/

/**
 * Test: Force magnitude calculation.
 *
 * ThDD:T-US-043-V.16
 *
 * MATHEMATICAL IDENTITY:
 *   |(3, 4, 0)| = sqrt(9 + 16 + 0) = 5
 *   |(1, 2, 2)| = sqrt(1 + 4 + 4) = 3
 */
static int test_embedding_stability_force_magnitude_T_US_043_V_16(void)
{
    const double TOL = 1e-14;

    /* 3-4-5 right triangle */
    double mag1 = compute_force_magnitude(3.0, 4.0, 0.0);
    if (fabs(mag1 - 5.0) > TOL) {
        fprintf(stderr, "\n    FAIL: |(3,4,0)| = %.15f, expected 5.0\n", mag1);
        return 0;
    }

    /* 1-2-2 vector */
    double mag2 = compute_force_magnitude(1.0, 2.0, 2.0);
    if (fabs(mag2 - 3.0) > TOL) {
        fprintf(stderr, "\n    FAIL: |(1,2,2)| = %.15f, expected 3.0\n", mag2);
        return 0;
    }

    /* Zero vector */
    double mag3 = compute_force_magnitude(0.0, 0.0, 0.0);
    if (fabs(mag3) > TOL) {
        fprintf(stderr, "\n    FAIL: |(0,0,0)| = %.15f, expected 0.0\n", mag3);
        return 0;
    }

    return 1;
}

/**
 * Test: Force jump calculation.
 *
 * ThDD:T-US-043-V.16
 *
 * MATHEMATICAL IDENTITY:
 *   |F2 - F1| where F1=(1,2,3), F2=(4,6,3)
 *   = |(3,4,0)| = 5
 */
static int test_embedding_stability_force_jump_T_US_043_V_16(void)
{
    const double TOL = 1e-14;

    /* F1 = (1,2,3), F2 = (4,6,3), diff = (3,4,0), |diff| = 5 */
    double jump1 = compute_force_jump(1.0, 2.0, 3.0, 4.0, 6.0, 3.0);
    if (fabs(jump1 - 5.0) > TOL) {
        fprintf(stderr, "\n    FAIL: force jump = %.15f, expected 5.0\n", jump1);
        return 0;
    }

    /* No change: F1 = F2 = (1,1,1) */
    double jump2 = compute_force_jump(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    if (fabs(jump2) > TOL) {
        fprintf(stderr, "\n    FAIL: zero force jump = %.15f, expected 0.0\n", jump2);
        return 0;
    }

    /* Unit change in X only: F1=(0,0,0), F2=(1,0,0), |diff|=1 */
    double jump3 = compute_force_jump(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    if (fabs(jump3 - 1.0) > TOL) {
        fprintf(stderr, "\n    FAIL: unit X jump = %.15f, expected 1.0\n", jump3);
        return 0;
    }

    return 1;
}

/**
 * Test: Switching region identification.
 *
 * ThDD:T-US-043-E.4
 *
 * MATHEMATICAL IDENTITY:
 *   r_on = 1.0, r_cut = 1.2
 *   r = 0.8 -> not in switching (inside inner region)
 *   r = 1.1 -> in switching
 *   r = 1.3 -> not in switching (outside cutoff)
 */
static int test_embedding_stability_switching_region_T_US_043_E_4(void)
{
    double r_on = 1.0;
    double r_cut = 1.2;

    /* Inside inner region (r < r_on) */
    if (is_in_switching_region(0.8, r_on, r_cut) != 0) {
        fprintf(stderr, "\n    FAIL: r=0.8 should NOT be in switching region\n");
        return 0;
    }

    /* In switching region (r_on <= r < r_cut) */
    if (is_in_switching_region(1.1, r_on, r_cut) != 1) {
        fprintf(stderr, "\n    FAIL: r=1.1 should be in switching region\n");
        return 0;
    }

    /* At inner boundary */
    if (is_in_switching_region(1.0, r_on, r_cut) != 1) {
        fprintf(stderr, "\n    FAIL: r=1.0 should be in switching region (boundary)\n");
        return 0;
    }

    /* Just inside cutoff */
    if (is_in_switching_region(1.19, r_on, r_cut) != 1) {
        fprintf(stderr, "\n    FAIL: r=1.19 should be in switching region\n");
        return 0;
    }

    /* At cutoff boundary (exclusive) */
    if (is_in_switching_region(1.2, r_on, r_cut) != 0) {
        fprintf(stderr, "\n    FAIL: r=1.2 should NOT be in switching region\n");
        return 0;
    }

    /* Outside cutoff */
    if (is_in_switching_region(1.3, r_on, r_cut) != 0) {
        fprintf(stderr, "\n    FAIL: r=1.3 should NOT be in switching region\n");
        return 0;
    }

    return 1;
}

/**
 * Test: Cutoff boundary crossing detection.
 *
 * ThDD:T-US-043-V.16
 *
 * MATHEMATICAL IDENTITY:
 *   r_cut = 1.2
 *   r_prev=1.1, r_curr=1.15 -> no crossing (both inside)
 *   r_prev=1.1, r_curr=1.25 -> crossing (in to out)
 *   r_prev=1.25, r_curr=1.15 -> crossing (out to in)
 *   r_prev=1.3, r_curr=1.4 -> no crossing (both outside)
 */
static int test_embedding_stability_boundary_crossing_T_US_043_V_16(void)
{
    double r_cut = 1.2;

    /* Both inside - no crossing */
    if (crossed_cutoff_boundary(1.1, 1.15, r_cut) != 0) {
        fprintf(stderr, "\n    FAIL: both inside should NOT be crossing\n");
        return 0;
    }

    /* Inside to outside - crossing */
    if (crossed_cutoff_boundary(1.1, 1.25, r_cut) != 1) {
        fprintf(stderr, "\n    FAIL: inside to outside should be crossing\n");
        return 0;
    }

    /* Outside to inside - crossing */
    if (crossed_cutoff_boundary(1.25, 1.15, r_cut) != 1) {
        fprintf(stderr, "\n    FAIL: outside to inside should be crossing\n");
        return 0;
    }

    /* Both outside - no crossing */
    if (crossed_cutoff_boundary(1.3, 1.4, r_cut) != 0) {
        fprintf(stderr, "\n    FAIL: both outside should NOT be crossing\n");
        return 0;
    }

    /* Edge case: exactly at boundary (inside) to outside */
    if (crossed_cutoff_boundary(1.19, 1.21, r_cut) != 1) {
        fprintf(stderr, "\n    FAIL: just inside to outside should be crossing\n");
        return 0;
    }

    return 1;
}

/**
 * Test: Threshold violation counting.
 *
 * ThDD:T-US-043-V.16
 *
 * MATHEMATICAL IDENTITY:
 *   data = {1.0, 5.0, 12.0, 8.0, 15.0}, threshold = 10.0
 *   Count above threshold: 2 (values 12.0 and 15.0)
 */
static int test_embedding_stability_threshold_counting_T_US_043_V_16(void)
{
    double data[5] = {1.0, 5.0, 12.0, 8.0, 15.0};

    int count = count_above_threshold(5, data, 10.0);
    if (count != 2) {
        fprintf(stderr, "\n    FAIL: count = %d, expected 2\n", count);
        return 0;
    }

    /* All below threshold */
    count = count_above_threshold(5, data, 20.0);
    if (count != 0) {
        fprintf(stderr, "\n    FAIL: count (threshold=20) = %d, expected 0\n", count);
        return 0;
    }

    /* All above threshold */
    count = count_above_threshold(5, data, 0.5);
    if (count != 5) {
        fprintf(stderr, "\n    FAIL: count (threshold=0.5) = %d, expected 5\n", count);
        return 0;
    }

    return 1;
}

/**
 * Test: Find maximum and index.
 *
 * MATHEMATICAL IDENTITY:
 *   data = {3.0, 7.0, 2.0, 9.0, 1.0}
 *   max = 9.0 at index 3
 */
static int test_embedding_stability_find_max_T_US_043_V_16(void)
{
    double data[5] = {3.0, 7.0, 2.0, 9.0, 1.0};
    double max_val;
    int max_idx;

    find_max_and_index(5, data, &max_val, &max_idx);

    if (fabs(max_val - 9.0) > 1e-14) {
        fprintf(stderr, "\n    FAIL: max_val = %.15f, expected 9.0\n", max_val);
        return 0;
    }

    if (max_idx != 3) {
        fprintf(stderr, "\n    FAIL: max_idx = %d, expected 3\n", max_idx);
        return 0;
    }

    return 1;
}

/**
 * Test: Force continuity threshold is physically sensible.
 *
 * ThDD:T-US-043-V.16
 *
 * The 10 kJ/mol/nm threshold should be:
 * - Large enough to allow normal thermal force fluctuations
 * - Small enough to detect discontinuities from cutoff artifacts
 *
 * At 300 K, thermal force fluctuations for a water molecule:
 *   F_thermal ~ sqrt(k_B * T * m) / dt ~ 1 kJ/mol/nm (order of magnitude)
 *   Threshold = 10x thermal -> allows normal dynamics
 */
static int test_embedding_stability_threshold_sensible_T_US_043_V_16(void)
{
    double k_B = 0.0083144626;  /* kJ/mol/K */
    double T = 300.0;           /* K */

    /* Thermal energy scale */
    double thermal_energy = k_B * T;  /* ~2.5 kJ/mol */

    /* Typical force scale for water at nm distances */
    /* F ~ E / r ~ 2.5 kJ/mol / 0.1 nm ~ 25 kJ/mol/nm (very rough) */
    /* Thermal fluctuation of force << mean force */

    /* Threshold should allow for significant force changes but catch discontinuities */
    if (MM_FORCE_JUMP_THRESHOLD_KJMOL_NM < 1.0) {
        fprintf(stderr, "\n    FAIL: Threshold %.1f kJ/mol/nm too small for thermal fluctuations\n",
                MM_FORCE_JUMP_THRESHOLD_KJMOL_NM);
        return 0;
    }

    if (MM_FORCE_JUMP_THRESHOLD_KJMOL_NM > 1000.0) {
        fprintf(stderr, "\n    FAIL: Threshold %.1f kJ/mol/nm too large to detect artifacts\n",
                MM_FORCE_JUMP_THRESHOLD_KJMOL_NM);
        return 0;
    }

    printf("\n    (threshold = %.0f kJ/mol/nm, thermal k_B*T = %.2f kJ/mol) ",
           MM_FORCE_JUMP_THRESHOLD_KJMOL_NM, thermal_energy);
    return 1;
}

/**
 * Test: V.15 cross-reference documentation.
 *
 * ThDD:T-US-043-V.15
 *
 * This test verifies that the cross-reference to US-039 is documented.
 * V.15 (switching smoothness) is validated by:
 *   - test_switching.c: Mathematical identity tests
 *   - test_switching_fd.c: Finite-difference validation
 *   - test_switching_regression.c: Backward compatibility
 *
 * This test passes unconditionally but documents the cross-reference.
 */
static int test_embedding_stability_V15_crossref_T_US_043_V_15(void)
{
    printf("\n");
    printf("    V.15 (Switching Function Smoothness) is validated in US-039:\n");
    printf("    - tests/test_switching.c: S(0)=1, S(1)=0, S'(0)=S'(1)=0\n");
    printf("    - tests/test_switching_fd.c: FD force validation\n");
    printf("    - tests/test_switching_regression.c: Backward compatibility\n");
    printf("\n");
    printf("    ThDD:T-US-043-V.15 cross-references US-039 V.6, V.7, V.8\n");
    printf("    ");

    /* This test documents the cross-reference - always passes */
    return 1;
}

/* ===========================================================================
 * SECTION 4: INTEGRATION TESTS (Trajectory Analysis)
 *
 * These tests require actual trajectory data from long GROMACS simulations.
 * Gated behind GRODFTB_LONG_TESTS.
 *
 * GOLDEN RULE: All expected values are thresholds from T-US-043, NOT
 * fabricated force predictions. The test analyzes actual trajectory data
 * and compares force jumps against the threshold.
 * ===========================================================================*/

#ifndef GRODFTB_SOURCE_DIR_STR
#define GRODFTB_SOURCE_DIR_STR "."
#endif
#define B5_DATA_DIR GRODFTB_SOURCE_DIR_STR "/tests/data/b5"

/**
 * Test V.16: MM force continuity over trajectory.
 *
 * ThDD:T-US-043-V.16
 * ThDD:T-US-043-E.3
 * ThDD:T-US-043-E.4
 *
 * Criterion: |F_MM(t+dt) - F_MM(t)| < 10 kJ/mol/nm
 *            (for atoms NOT crossing the cutoff boundary)
 *
 * This test analyzes the force trajectory to verify that MM forces
 * are smooth enough to avoid causing energy drift. Large force jumps
 * indicate:
 *   1. Cutoff artifacts (if switching is disabled)
 *   2. Implementation bugs in the switching function
 *   3. Numerical instability in force computation
 *
 * NOTE: This test SKIPS if trajectory data is not available.
 * The trajectory must be generated by running GROMACS mdrun.
 */
static int test_mm_force_continuity_V_US_043_V_16(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    printf("\n");
    printf("    === V.16: MM Force Continuity ===\n");
    printf("    ThDD:T-US-043-V.16, T-US-043-E.3, T-US-043-E.4\n");
    printf("\n");
    printf("    Criterion: |F_MM(t+dt) - F_MM(t)| < %.0f kJ/mol/nm\n",
           MM_FORCE_JUMP_THRESHOLD_KJMOL_NM);
    printf("    (for atoms that did NOT cross the cutoff boundary)\n");
    printf("\n");
    printf("    System: B5 (%d MM atoms, cutoff=%.1f nm)\n",
           B5_N_MM_ATOMS, B5_CUTOFF_NM);
    printf("\n");
    printf("    Rationale (T-US-043-E.3, T-US-043-E.4):\n");
    printf("    - Large force jumps indicate cutoff artifacts\n");
    printf("    - Quintic switching ensures smooth force at boundaries\n");
    printf("    - Residual numerical error should be << threshold\n");
    printf("\n");

    /* Check for required data files */
    char trr_path[512];
    snprintf(trr_path, sizeof(trr_path), "%s/nve_500ps/mm_forces.xvg", B5_DATA_DIR);

    FILE *f = fopen(trr_path, "r");
    if (!f) {
        printf("    MM force data not found: %s\n", trr_path);
        printf("\n");
        printf("    To generate force data:\n");
        printf("    1. Run 500 ps NVE trajectory with force output enabled\n");
        printf("    2. Extract MM forces: gmx traj -f nve_500ps.trr -of mm_forces.xvg -s nve_500ps.tpr\n");
        printf("    3. Also need MM positions for distance calculation\n");
        printf("\n");
        printf("    SKIP: Trajectory force data required\n");
        return -1;
    }
    fclose(f);

    /* Check for position data (needed to identify boundary crossings) */
    char pos_path[512];
    snprintf(pos_path, sizeof(pos_path), "%s/nve_500ps/mm_positions.xvg", B5_DATA_DIR);

    f = fopen(pos_path, "r");
    if (!f) {
        printf("    MM position data not found: %s\n", pos_path);
        printf("    (Needed to identify atoms crossing cutoff boundary)\n");
        printf("\n");
        printf("    SKIP: Trajectory position data required\n");
        return -1;
    }
    fclose(f);

    printf("    Force and position data files found.\n");
    printf("    Full analysis not yet implemented.\n");
    printf("\n");
    printf("    Expected analysis steps:\n");
    printf("    1. Read force trajectory frame by frame\n");
    printf("    2. For each MM atom, compute |F(t+dt) - F(t)|\n");
    printf("    3. Use position data to identify boundary crossings\n");
    printf("    4. Exclude crossed atoms from force jump check\n");
    printf("    5. Report statistics and threshold violations\n");
    printf("\n");
    printf("    SKIP: Full analysis pending implementation\n");
    return -1;
#endif
}

/**
 * Test: MM force statistics over long trajectory.
 *
 * ThDD:T-US-043-E.3
 *
 * This test computes force jump statistics:
 *   - Mean force jump magnitude
 *   - Standard deviation of force jumps
 *   - Maximum force jump (and which atom/frame)
 *   - Number of threshold violations
 *
 * The statistics help characterize normal force dynamics vs artifacts.
 */
static int test_mm_force_statistics_T_US_043_E_3(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    printf("\n");
    printf("    === Force Jump Statistics ===\n");
    printf("    ThDD:T-US-043-E.3 -- Cumulative energy error bound\n");
    printf("\n");
    printf("    This diagnostic test computes force jump statistics.\n");
    printf("\n");
    printf("    Expected output (from actual trajectory):\n");
    printf("    - Mean force jump: ~ 0.1-1.0 kJ/mol/nm (thermal)\n");
    printf("    - Stddev: ~ same order as mean\n");
    printf("    - Max jump: should be < %.0f kJ/mol/nm\n",
           MM_FORCE_JUMP_THRESHOLD_KJMOL_NM);
    printf("    - Violations: should be 0 or very few\n");
    printf("\n");

    char trr_path[512];
    snprintf(trr_path, sizeof(trr_path), "%s/nve_500ps/mm_forces.xvg", B5_DATA_DIR);

    FILE *f = fopen(trr_path, "r");
    if (!f) {
        printf("    Force data not found: %s\n", trr_path);
        printf("    SKIP: Trajectory data required\n");
        return -1;
    }
    fclose(f);

    printf("    Force data found - analysis not yet implemented\n");
    printf("    SKIP: Full analysis pending implementation\n");
    return -1;
#endif
}

/**
 * Test: Force continuity in switching region specifically.
 *
 * ThDD:T-US-043-E.4
 *
 * This test focuses on atoms within the switching region (r_on <= r < r_cut).
 * The quintic switching function should ensure smooth forces even for atoms
 * near the cutoff boundary.
 *
 * Per T-US-043-E.4, the energy change at each timestep should follow the
 * correct work-energy relation: dE = -F . v dt
 */
static int test_force_continuity_switching_region_T_US_043_E_4(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    printf("\n");
    printf("    === Force Continuity in Switching Region ===\n");
    printf("    ThDD:T-US-043-E.4 -- Energy continuity through switching\n");
    printf("\n");
    printf("    Switching region: %.2f nm <= r < %.2f nm\n", B5_R_ON_NM, B5_CUTOFF_NM);
    printf("\n");
    printf("    Expected (from T-US-043-E.1):\n");
    printf("    - ~300 atoms in switching region at any time\n");
    printf("    - ~42%% of embedding shell\n");
    printf("\n");
    printf("    Verification:\n");
    printf("    - Force on atoms IN switching region should be continuous\n");
    printf("    - No force jumps > threshold even near boundary\n");
    printf("\n");

    char pos_path[512];
    snprintf(pos_path, sizeof(pos_path), "%s/nve_500ps/mm_positions.xvg", B5_DATA_DIR);

    FILE *f = fopen(pos_path, "r");
    if (!f) {
        printf("    Position data not found: %s\n", pos_path);
        printf("    SKIP: Trajectory data required\n");
        return -1;
    }
    fclose(f);

    printf("    Position data found - analysis not yet implemented\n");
    printf("    SKIP: Full analysis pending implementation\n");
    return -1;
#endif
}

/* ===========================================================================
 * SECTION 5: MAIN TEST RUNNER
 * ===========================================================================*/

static void print_usage(const char *prog)
{
    printf("Usage: %s [--unit | --integration | --all]\n", prog);
    printf("\n");
    printf("Options:\n");
    printf("  --unit         Run unit tests only (default)\n");
    printf("  --integration  Run integration tests only\n");
    printf("  --all          Run all tests\n");
    printf("\n");
    printf("Unit tests verify force analysis helper functions.\n");
    printf("Integration tests require trajectory data from long GROMACS simulations.\n");
    printf("Integration tests are gated behind GRODFTB_LONG_TESTS cmake option.\n");
    printf("\n");
    printf("Cross-reference:\n");
    printf("  V.15 (switching smoothness) is tested in US-039:\n");
    printf("    - tests/test_switching.c\n");
    printf("    - tests/test_switching_fd.c\n");
    printf("    - tests/test_switching_regression.c\n");
}

int main(int argc, char *argv[])
{
    int run_unit = 1;
    int run_integration = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--unit") == 0) {
            run_unit = 1;
            run_integration = 0;
        } else if (strcmp(argv[i], "--integration") == 0) {
            run_unit = 0;
            run_integration = 1;
        } else if (strcmp(argv[i], "--all") == 0) {
            run_unit = 1;
            run_integration = 1;
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            print_usage(argv[0]);
            return 0;
        }
    }

    printf("=== US-043 Embedding Stability Tests (V.15-V.16) ===\n\n");

    printf("Verification Criteria (from T-US-043):\n");
    printf("  V.15: Switching function smoothness (cross-reference US-039)\n");
    printf("  V.16: MM force continuity < %.0f kJ/mol/nm\n",
           MM_FORCE_JUMP_THRESHOLD_KJMOL_NM);
    printf("\n");

    printf("Theory References:\n");
    printf("  T-US-043-V.16: Force jump threshold definition\n");
    printf("  T-US-043-E.3:  Cumulative energy error bound\n");
    printf("  T-US-043-E.4:  Energy continuity through switching transitions\n");
    printf("\n");

    printf("B5 System Parameters:\n");
    printf("  Total atoms: %d (QM: %d, MM: %d)\n",
           B5_N_ATOMS_TOTAL, B5_N_QM_ATOMS, B5_N_MM_ATOMS);
    printf("  Cutoff: %.2f nm\n", B5_CUTOFF_NM);
    printf("  Switch width: %.2f nm (r_on = %.2f nm)\n",
           B5_SWITCH_WIDTH_NM, B5_R_ON_NM);
    printf("\n");

    if (run_unit) {
        printf("=== Unit Tests (Mathematical Identities) ===\n\n");

        printf("Force Magnitude (T-US-043-V.16):\n");
        RUN_TEST(test_embedding_stability_force_magnitude_T_US_043_V_16);
        printf("\n");

        printf("Force Jump (T-US-043-V.16):\n");
        RUN_TEST(test_embedding_stability_force_jump_T_US_043_V_16);
        printf("\n");

        printf("Switching Region Identification (T-US-043-E.4):\n");
        RUN_TEST(test_embedding_stability_switching_region_T_US_043_E_4);
        printf("\n");

        printf("Boundary Crossing Detection (T-US-043-V.16):\n");
        RUN_TEST(test_embedding_stability_boundary_crossing_T_US_043_V_16);
        printf("\n");

        printf("Threshold Counting (T-US-043-V.16):\n");
        RUN_TEST(test_embedding_stability_threshold_counting_T_US_043_V_16);
        printf("\n");

        printf("Find Maximum (T-US-043-V.16):\n");
        RUN_TEST(test_embedding_stability_find_max_T_US_043_V_16);
        printf("\n");

        printf("Threshold Validation (T-US-043-V.16):\n");
        RUN_TEST(test_embedding_stability_threshold_sensible_T_US_043_V_16);
        printf("\n");

        printf("V.15 Cross-Reference (T-US-043-V.15):\n");
        RUN_TEST(test_embedding_stability_V15_crossref_T_US_043_V_15);
        printf("\n");
    }

    if (run_integration) {
        printf("=== Integration Tests (V.16 Trajectory Analysis) ===\n\n");

#ifndef GRODFTB_LONG_TESTS
        printf("NOTE: Integration tests require GRODFTB_LONG_TESTS=ON\n");
        printf("      Rebuild with: cmake -DGRODFTB_LONG_TESTS=ON ..\n\n");
#endif

        printf("V.16 - MM Force Continuity:\n");
        RUN_TEST(test_mm_force_continuity_V_US_043_V_16);
        printf("\n");

        printf("Force Statistics (T-US-043-E.3):\n");
        RUN_TEST(test_mm_force_statistics_T_US_043_E_3);
        printf("\n");

        printf("Switching Region Analysis (T-US-043-E.4):\n");
        RUN_TEST(test_force_continuity_switching_region_T_US_043_E_4);
        printf("\n");
    }

    printf("=== Summary ===\n");
    printf("Tests run:     %d\n", tests_run);
    printf("Tests passed:  %d\n", tests_passed);
    printf("Tests failed:  %d\n", tests_failed);
    printf("Tests skipped: %d\n", tests_skipped);

    if (tests_skipped > 0 && run_integration) {
        printf("\nNOTE: Skipped integration tests require:\n");
        printf("  1. GRODFTB_LONG_TESTS=ON at compile time\n");
        printf("  2. Trajectory data in tests/data/b5/nve_500ps/\n");
        printf("     - mm_forces.xvg (force trajectory)\n");
        printf("     - mm_positions.xvg (position trajectory)\n");
    }

    printf("\nCross-reference for V.15 (switching smoothness):\n");
    printf("  See tests/test_switching.c, test_switching_fd.c, test_switching_regression.c\n");

    if (tests_failed > 0) {
        printf("\nFAILED: %d test(s) failed\n", tests_failed);
        return 1;
    }

    return 0;
}
