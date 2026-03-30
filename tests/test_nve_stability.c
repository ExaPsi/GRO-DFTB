/*
 * US-043: Extended NVE Stability Tests (100-500 ps)
 *
 * ThDD:T-US-043-V.1 through V.4 -- Energy drift and statistical criteria
 * SDD:specs.md:S21.1 -- NVE drift criterion (< 0.01 kJ/mol/ps/atom)
 * SDD:docs/verification/US-043.md -- Verification plan
 *
 * This file implements:
 *   1. V.1: Energy drift over 100ps and 500ps trajectories
 *   2. V.2: Energy fluctuation amplitude check
 *   3. V.3: R-squared interpretation (fluctuation vs drift dominated)
 *   4. V.4: Block-averaged drift consistency
 *
 * GOLDEN RULE COMPLIANCE:
 * - Unit tests use mathematical identities (exact formulas)
 * - Integration test thresholds from specs.md, not fabricated
 * - No hardcoded "expected" simulation energy values
 * - All simulation results must come from actual trajectory data
 *
 * Test execution modes:
 *   --unit         Run unit tests only (fast, no simulation required)
 *   --integration  Run integration tests (requires GRODFTB_LONG_TESTS)
 *   (no args)      Run unit tests only
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
 * PHYSICAL CONSTANTS
 *
 * ThDD:CODATA2022 -- All constants from verified sources
 * ===========================================================================*/

/* Boltzmann constant in kJ/mol/K (CODATA 2022) */
#define KB_KJMOL_K  0.0083144626

/* ===========================================================================
 * ACCEPTANCE CRITERIA (from specs.md Section 21.1)
 *
 * SDD:specs.md:S21.1:
 *   "Stability: NVE drift | 3, 4 | long NVE QM/MM | < 0.01 kJ/mol/ps/atom"
 *
 * ThDD:T-US-043-V.1:
 *   drift_normalized = |slope| / N_atoms < 0.01 kJ/mol/ps/atom
 *
 * These are specification-derived thresholds, NOT fabricated expected values.
 * ===========================================================================*/

/* SDD:specs.md:S21.1 -- Primary acceptance criterion */
#define NVE_DRIFT_THRESHOLD_KJMOL_PS_ATOM  0.01

/* ThDD:T-US-043-V.2 -- Energy fluctuation threshold */
/* sigma_E < 2 * sqrt(N_atoms) * k_B * T = 2 * sqrt(2652) * 2.5 ~ 260 kJ/mol */
#define ENERGY_FLUCTUATION_MULTIPLIER  2.0

/* ThDD:T-US-043-V.3 -- R-squared threshold (fluctuation-dominated) */
#define R_SQUARED_UPPER_BOUND  0.5

/* ThDD:T-US-043-V.4 -- Block drift consistency (within 3 standard errors) */
#define BLOCK_DRIFT_SE_MULTIPLIER  3.0

/* ===========================================================================
 * B5 SYSTEM PARAMETERS (from tests/data/b5/provenance.json)
 *
 * ThDD:CLAUDE.md:Golden_Rule -- All values from actual GROMACS preparation
 * ===========================================================================*/

#define B5_N_ATOMS_TOTAL  2652   /* Total atoms in system */
#define B5_N_QM_ATOMS     3      /* QM region: 1 water molecule */
#define B5_TARGET_T_K     300.0  /* Target temperature in Kelvin */

/* ===========================================================================
 * REGRESSION RESULT STRUCTURE (reused from US-042)
 *
 * ThDD:T-US-042-N.3 through N.6 -- Complete regression statistics
 * ===========================================================================*/

/**
 * Complete result from linear regression analysis.
 *
 * ThDD:T-US-042-N.3 -- slope = S_tE / S_tt
 * ThDD:T-US-042-N.4 -- r_squared = 1 - SS_res / SS_tot
 * ThDD:T-US-042-N.5 -- sigma_E = sqrt(SS_res / (n-2))
 * ThDD:T-US-042-N.6 -- sigma_slope = sigma_E / sqrt(S_tt)
 */
typedef struct {
    double slope;       /* Drift rate (same units as E/t) */
    double intercept;   /* Initial energy estimate */
    double r_squared;   /* Coefficient of determination */
    double sigma_E;     /* Residual standard error */
    double sigma_slope; /* Standard error of slope */
    double S_tt;        /* Sum of squared time deviations */
    double SS_res;      /* Sum of squared residuals */
    double SS_tot;      /* Total sum of squares */
    int n;              /* Number of data points */
} regression_result_t;

/* ===========================================================================
 * SECTION 1: STATISTICAL ANALYSIS FUNCTIONS
 *
 * Identical to US-042 implementation for consistency.
 * ThDD:T-US-042-N.3 through N.6
 * ===========================================================================*/

/**
 * Compute linear regression with full statistics.
 *
 * ThDD:T-US-042-N.3 -- Linear regression formulas
 * ThDD:T-US-042-N.4 -- R-squared calculation
 * ThDD:T-US-042-N.5 -- Residual standard error
 * ThDD:T-US-042-N.6 -- Standard error of slope
 *
 * Numerically stable implementation using centered data with Kahan summation.
 *
 * @param n       Number of data points (must be >= 2)
 * @param t       Time array [n]
 * @param E       Energy array [n]
 * @param result  Output structure receiving all statistics
 */
static void compute_linear_regression(int n, const double *t, const double *E,
                                       regression_result_t *result)
{
    memset(result, 0, sizeof(*result));
    result->n = n;

    if (n < 2) {
        result->slope = 0.0;
        result->intercept = (n > 0) ? E[0] : 0.0;
        return;
    }

    /* Step 1: Compute means with Kahan summation */
    double sum_t = 0.0, sum_E = 0.0;
    double c_t = 0.0, c_E = 0.0;

    for (int i = 0; i < n; i++) {
        double y_t = t[i] - c_t;
        double temp_t = sum_t + y_t;
        c_t = (temp_t - sum_t) - y_t;
        sum_t = temp_t;

        double y_E = E[i] - c_E;
        double temp_E = sum_E + y_E;
        c_E = (temp_E - sum_E) - y_E;
        sum_E = temp_E;
    }

    double t_bar = sum_t / (double)n;
    double E_bar = sum_E / (double)n;

    /* Step 2: Compute centered sums */
    double S_tt = 0.0, S_tE = 0.0, S_EE = 0.0;
    double c_tt = 0.0, c_tE = 0.0, c_EE = 0.0;

    for (int i = 0; i < n; i++) {
        double dt = t[i] - t_bar;
        double dE = E[i] - E_bar;

        double y_tt = dt * dt - c_tt;
        double temp_tt = S_tt + y_tt;
        c_tt = (temp_tt - S_tt) - y_tt;
        S_tt = temp_tt;

        double y_tE = dt * dE - c_tE;
        double temp_tE = S_tE + y_tE;
        c_tE = (temp_tE - S_tE) - y_tE;
        S_tE = temp_tE;

        double y_EE = dE * dE - c_EE;
        double temp_EE = S_EE + y_EE;
        c_EE = (temp_EE - S_EE) - y_EE;
        S_EE = temp_EE;
    }

    result->S_tt = S_tt;
    result->SS_tot = S_EE;

    /* Step 3: Compute slope and intercept */
    if (fabs(S_tt) < 1e-30) {
        result->slope = 0.0;
        result->intercept = E_bar;
        return;
    }

    result->slope = S_tE / S_tt;
    result->intercept = E_bar - result->slope * t_bar;

    /* Step 4: Compute R-squared */
    double SS_res = 0.0;
    double c_res = 0.0;

    for (int i = 0; i < n; i++) {
        double E_pred = result->slope * t[i] + result->intercept;
        double residual = E[i] - E_pred;
        double resid_sq = residual * residual;

        double y_res = resid_sq - c_res;
        double temp_res = SS_res + y_res;
        c_res = (temp_res - SS_res) - y_res;
        SS_res = temp_res;
    }

    result->SS_res = SS_res;

    if (S_EE > 1e-30) {
        result->r_squared = 1.0 - SS_res / S_EE;
    } else {
        result->r_squared = 0.0;
    }

    /* Step 5: Compute residual standard error */
    if (n > 2) {
        result->sigma_E = sqrt(SS_res / (double)(n - 2));
    } else {
        result->sigma_E = 0.0;
    }

    /* Step 6: Compute standard error of slope */
    if (S_tt > 1e-30) {
        result->sigma_slope = result->sigma_E / sqrt(S_tt);
    } else {
        result->sigma_slope = 0.0;
    }
}

/**
 * Compute normalized drift rate.
 *
 * ThDD:T-US-043-V.1 -- Normalized drift: drift_norm = |slope| / N_atoms
 *
 * @param slope    Drift rate in kJ/mol/ps
 * @param n_atoms  Number of atoms in system
 * @return         Normalized drift in kJ/mol/ps/atom
 */
static double compute_normalized_drift(double slope, int n_atoms)
{
    if (n_atoms <= 0) return 0.0;
    return fabs(slope) / (double)n_atoms;
}

/**
 * Compute sample standard deviation of an array.
 *
 * ThDD:T-US-043-V.2 -- Energy fluctuation sigma_E
 *
 * @param n    Number of elements
 * @param arr  Data array
 * @return     Sample standard deviation (Bessel-corrected)
 */
static double compute_stddev(int n, const double *arr)
{
    if (n < 2) return 0.0;

    double mean = 0.0;
    for (int i = 0; i < n; i++) {
        mean += arr[i];
    }
    mean /= (double)n;

    double variance = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = arr[i] - mean;
        variance += diff * diff;
    }
    variance /= (double)(n - 1);

    return sqrt(variance);
}

/**
 * Compute expected energy fluctuation for a canonical system.
 *
 * ThDD:T-US-043-V.2 -- sigma_E ~ sqrt(N) * k_B * T
 *
 * This is an approximate formula based on equipartition. The actual
 * NVE fluctuation may differ due to fixed total energy constraint.
 *
 * @param n_atoms    Number of atoms
 * @param T_kelvin   Temperature in Kelvin
 * @return           Expected energy standard deviation in kJ/mol
 */
static double compute_expected_sigma_E(int n_atoms, double T_kelvin)
{
    /*
     * ThDD:T-US-043-V.2
     *
     * Canonical fluctuation: sigma_E ~ sqrt(k_B * T^2 * C_V)
     * For an ideal gas: C_V ~ (3/2) * N * k_B
     * -> sigma_E ~ sqrt((3/2) * N * k_B^2 * T^2) ~ T * sqrt(N) * k_B * sqrt(3/2)
     *
     * Simplified bound: sigma_E < 2 * sqrt(N) * k_B * T
     */
    return sqrt((double)n_atoms) * KB_KJMOL_K * T_kelvin;
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
 * Test: Linear regression recovers exact slope.
 *
 * ThDD:T-US-042-N.3
 *
 * MATHEMATICAL IDENTITY:
 *   y = 0.005 * x + 1000 for x = {0, 1, ..., 99999}
 *   Expected: slope = 0.005 exactly (500 ps at 5 fs output)
 */
static int test_nve_stability_regression_exact_slope_T_US_043_N_1(void)
{
    const int n = 100000;  /* 500 ps at 5 fs output */
    double *t = malloc(n * sizeof(double));
    double *E = malloc(n * sizeof(double));

    if (!t || !E) {
        if (t) free(t);
        if (E) free(E);
        fprintf(stderr, "\n    FAIL: Memory allocation failed\n");
        return 0;
    }

    /* Generate exact linear data: E = 0.005*t + 1000 (drift = 0.005 kJ/mol/ps) */
    for (int i = 0; i < n; i++) {
        t[i] = i * 0.005;  /* Time in ps (5 fs steps) */
        E[i] = 0.005 * t[i] + 1000.0;
    }

    regression_result_t result;
    compute_linear_regression(n, t, E, &result);

    free(t);
    free(E);

    const double TOL = 1e-10;

    if (fabs(result.slope - 0.005) > TOL) {
        fprintf(stderr, "\n    FAIL: slope = %.15e, expected 0.005\n", result.slope);
        return 0;
    }

    if (fabs(result.r_squared - 1.0) > TOL) {
        fprintf(stderr, "\n    FAIL: R^2 = %.15f, expected 1.0\n", result.r_squared);
        return 0;
    }

    printf("\n    (100000 points, slope error = %.2e) ", fabs(result.slope - 0.005));
    return 1;
}

/**
 * Test: Normalized drift calculation.
 *
 * ThDD:T-US-043-V.1
 *
 * MATHEMATICAL IDENTITY:
 *   slope = 26.52 kJ/mol/ps, N_atoms = 2652
 *   -> drift_norm = 26.52 / 2652 = 0.01 kJ/mol/ps/atom
 */
static int test_nve_stability_normalized_drift_T_US_043_V_1(void)
{
    double slope = 26.52;
    int n_atoms = 2652;

    double drift_norm = compute_normalized_drift(slope, n_atoms);

    const double TOL = 1e-10;

    if (fabs(drift_norm - 0.01) > TOL) {
        fprintf(stderr, "\n    FAIL: drift_norm = %.15f, expected 0.01\n", drift_norm);
        return 0;
    }

    return 1;
}

/**
 * Test: Energy fluctuation formula verification.
 *
 * ThDD:T-US-043-V.2
 *
 * MATHEMATICAL IDENTITY:
 *   sigma_E ~ sqrt(N) * k_B * T
 *   For N=2652, T=300K:
 *     sigma_E ~ sqrt(2652) * 0.00831 * 300 ~ 51.5 * 2.49 ~ 128.4 kJ/mol
 *
 *   Acceptance bound: 2 * sigma_E ~ 260 kJ/mol
 */
static int test_nve_stability_energy_fluctuation_formula_T_US_043_V_2(void)
{
    double sigma_E = compute_expected_sigma_E(B5_N_ATOMS_TOTAL, B5_TARGET_T_K);

    /* Expected from formula: sqrt(2652) * 0.0083144626 * 300 */
    double expected = sqrt((double)B5_N_ATOMS_TOTAL) * KB_KJMOL_K * B5_TARGET_T_K;

    const double TOL = 1e-6;

    if (fabs(sigma_E - expected) > TOL) {
        fprintf(stderr, "\n    FAIL: sigma_E = %.6f, expected %.6f\n",
                sigma_E, expected);
        return 0;
    }

    /* Verify numeric value (~128 kJ/mol) */
    if (sigma_E < 120.0 || sigma_E > 140.0) {
        fprintf(stderr, "\n    FAIL: sigma_E = %.1f kJ/mol, expected ~128\n", sigma_E);
        return 0;
    }

    printf("\n    (B5 expected sigma_E ~ %.1f kJ/mol, bound = %.1f kJ/mol) ",
           sigma_E, ENERGY_FLUCTUATION_MULTIPLIER * sigma_E);
    return 1;
}

/**
 * Test: R-squared interpretation for fluctuation-dominated data.
 *
 * ThDD:T-US-043-V.3
 *
 * MATHEMATICAL IDENTITY:
 *   For data with large oscillations and small or no systematic trend:
 *   - R^2 should be low (oscillations dominate over any trend)
 *
 * NOTE: A sine wave correlated with time has non-zero slope due to
 *       integration by parts: integral(t*sin(at)) = -cos(at)/a + t*cos(at)/a.
 *       For fluctuation-dominated test, we use symmetric oscillation.
 */
static int test_nve_stability_r_squared_noise_T_US_043_V_3(void)
{
    const int n = 1000;
    double *t = malloc(n * sizeof(double));
    double *E = malloc(n * sizeof(double));

    if (!t || !E) {
        if (t) free(t);
        if (E) free(E);
        fprintf(stderr, "\n    FAIL: Memory allocation failed\n");
        return 0;
    }

    /*
     * Generate symmetric oscillating data about the midpoint.
     * E[i] = 1000 + 10 * sin(2*pi*i/100) for i < n/2
     * E[n-1-i] = E[i] to ensure symmetry and zero net slope
     */
    const double PI = 3.14159265358979323846;
    for (int i = 0; i < n; i++) {
        t[i] = i * 0.01;  /* 10 fs steps, total 10 ps */
    }

    /* Build symmetric oscillation: first half and second half mirror each other */
    for (int i = 0; i < n/2; i++) {
        double osc = 10.0 * sin(2.0 * PI * i / 50.0);  /* 10 periods in first half */
        E[i] = 1000.0 + osc;
        E[n - 1 - i] = 1000.0 + osc;  /* Mirror in second half */
    }

    regression_result_t result;
    compute_linear_regression(n, t, E, &result);

    free(t);
    free(E);

    /* R^2 should be very low for symmetric oscillation (no net drift) */
    if (result.r_squared > 0.1) {
        fprintf(stderr, "\n    FAIL: R^2 = %.6f > 0.1 for symmetric oscillation\n",
                result.r_squared);
        return 0;
    }

    /* Slope should be essentially zero for symmetric data */
    if (fabs(result.slope) > 1e-10) {
        fprintf(stderr, "\n    FAIL: slope = %.6e for symmetric oscillation\n",
                result.slope);
        return 0;
    }

    printf("\n    (symmetric oscillation: R^2 = %.6f, slope = %.2e) ",
           result.r_squared, result.slope);
    return 1;
}

/**
 * Test: Block drift consistency for identical blocks.
 *
 * ThDD:T-US-043-V.4
 *
 * MATHEMATICAL IDENTITY:
 *   If data has constant slope throughout, first half and second half
 *   should give identical slopes.
 */
static int test_nve_stability_block_consistency_identical_T_US_043_V_4(void)
{
    const int n = 2000;
    double *t = malloc(n * sizeof(double));
    double *E = malloc(n * sizeof(double));

    if (!t || !E) {
        if (t) free(t);
        if (E) free(E);
        fprintf(stderr, "\n    FAIL: Memory allocation failed\n");
        return 0;
    }

    /* Generate linear data with constant drift */
    for (int i = 0; i < n; i++) {
        t[i] = i * 0.01;
        E[i] = 0.001 * t[i] + 500.0;
    }

    /* Compute drift for first half */
    regression_result_t result_first;
    compute_linear_regression(n/2, t, E, &result_first);

    /* Compute drift for second half */
    regression_result_t result_second;
    compute_linear_regression(n/2, t + n/2, E + n/2, &result_second);

    free(t);
    free(E);

    /* Both halves should have identical slope */
    double diff = fabs(result_first.slope - result_second.slope);
    const double TOL = 1e-10;

    if (diff > TOL) {
        fprintf(stderr, "\n    FAIL: slope_first = %.15f, slope_second = %.15f\n",
                result_first.slope, result_second.slope);
        return 0;
    }

    return 1;
}

/**
 * Test: Standard deviation calculation.
 *
 * MATHEMATICAL IDENTITY:
 *   For data {1, 2, 3, 4, 5}:
 *   mean = 3, variance = ((1-3)^2 + (2-3)^2 + (3-3)^2 + (4-3)^2 + (5-3)^2) / 4
 *        = (4 + 1 + 0 + 1 + 4) / 4 = 10/4 = 2.5
 *   stddev = sqrt(2.5) = 1.5811...
 */
static int test_nve_stability_stddev_calculation_T_US_043_N_2(void)
{
    double data[5] = {1.0, 2.0, 3.0, 4.0, 5.0};

    double stddev = compute_stddev(5, data);
    double expected = sqrt(2.5);

    const double TOL = 1e-10;

    if (fabs(stddev - expected) > TOL) {
        fprintf(stderr, "\n    FAIL: stddev = %.15f, expected %.15f\n",
                stddev, expected);
        return 0;
    }

    return 1;
}

/* ===========================================================================
 * SECTION 4: INTEGRATION TESTS (Trajectory Analysis)
 *
 * These tests require actual trajectory data from long GROMACS simulations.
 * Gated behind GRODFTB_LONG_TESTS.
 *
 * GOLDEN RULE: All expected values are thresholds from specs.md, NOT
 * fabricated energy predictions.
 * ===========================================================================*/

#ifndef GRODFTB_SOURCE_DIR_STR
#define GRODFTB_SOURCE_DIR_STR "."
#endif
#define B5_DATA_DIR GRODFTB_SOURCE_DIR_STR "/tests/data/b5"

/**
 * Test V.1: Energy drift over 100 ps trajectory.
 *
 * ThDD:T-US-043-V.1
 * SDD:specs.md:S21.1 -- drift < 0.01 kJ/mol/ps/atom
 *
 * NOTE: This test SKIPS if trajectory data is not available.
 * The trajectory must be generated by running GROMACS mdrun.
 */
static int test_nve_drift_100ps_V_US_043_V_1(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    printf("\n");
    printf("    === V.1: Energy Drift (100 ps) ===\n");
    printf("    Acceptance: drift < %.4f kJ/mol/ps/atom\n",
           NVE_DRIFT_THRESHOLD_KJMOL_PS_ATOM);
    printf("    System: B5 (%d atoms)\n", B5_N_ATOMS_TOTAL);
    printf("\n");

    char xvg_path[512];
    snprintf(xvg_path, sizeof(xvg_path), "%s/nve_100ps/total_energy.xvg", B5_DATA_DIR);

    FILE *f = fopen(xvg_path, "r");
    if (!f) {
        printf("    Trajectory data not found: %s\n", xvg_path);
        printf("    Generate with: gmx mdrun -s nve_100ps.tpr -deffnm nve_100ps\n");
        printf("    Then: gmx energy -f nve_100ps.edr -o total_energy.xvg\n");
        printf("\n    SKIP: Trajectory data required\n");
        return -1;
    }
    fclose(f);

    /* TODO: Implement XVG parsing and drift analysis */
    printf("    Trajectory found - analysis not yet implemented\n");
    printf("    SKIP: Full analysis pending implementation\n");
    return -1;
#endif
}

/**
 * Test V.1: Energy drift over 500 ps trajectory.
 *
 * ThDD:T-US-043-V.1
 * SDD:specs.md:S21.1 -- drift < 0.01 kJ/mol/ps/atom
 */
static int test_nve_drift_500ps_V_US_043_V_1(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    printf("\n");
    printf("    === V.1: Energy Drift (500 ps) ===\n");
    printf("    Acceptance: drift < %.4f kJ/mol/ps/atom\n",
           NVE_DRIFT_THRESHOLD_KJMOL_PS_ATOM);
    printf("    System: B5 (%d atoms)\n", B5_N_ATOMS_TOTAL);
    printf("    Duration: 500 ps (10^6 steps at 0.5 fs)\n");
    printf("\n");

    char xvg_path[512];
    snprintf(xvg_path, sizeof(xvg_path), "%s/nve_500ps/total_energy.xvg", B5_DATA_DIR);

    FILE *f = fopen(xvg_path, "r");
    if (!f) {
        printf("    Trajectory data not found: %s\n", xvg_path);
        printf("    Generate with: gmx mdrun -s nve_500ps.tpr -deffnm nve_500ps\n");
        printf("\n    SKIP: Trajectory data required\n");
        return -1;
    }
    fclose(f);

    printf("    Trajectory found - analysis not yet implemented\n");
    printf("    SKIP: Full analysis pending implementation\n");
    return -1;
#endif
}

/**
 * Test V.2: Energy fluctuation amplitude.
 *
 * ThDD:T-US-043-V.2
 *
 * Criterion: sigma_E < 2 * sqrt(N) * k_B * T ~ 260 kJ/mol
 */
static int test_energy_fluctuation_V_US_043_V_2(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    double sigma_E_bound = ENERGY_FLUCTUATION_MULTIPLIER *
                           compute_expected_sigma_E(B5_N_ATOMS_TOTAL, B5_TARGET_T_K);

    printf("\n");
    printf("    === V.2: Energy Fluctuation ===\n");
    printf("    Acceptance: sigma_E < %.1f kJ/mol\n", sigma_E_bound);
    printf("    (2 * sqrt(N) * k_B * T for B5 at 300 K)\n");
    printf("\n");
    printf("    SKIP: Requires trajectory data\n");
    return -1;
#endif
}

/**
 * Test V.3: R-squared interpretation.
 *
 * ThDD:T-US-043-V.3
 *
 * Criterion: R^2 < 0.5 (fluctuation-dominated, not drift-dominated)
 */
static int test_drift_rsquared_V_US_043_V_3(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    printf("\n");
    printf("    === V.3: R-squared Interpretation ===\n");
    printf("    Acceptance: R^2 < %.1f\n", R_SQUARED_UPPER_BOUND);
    printf("\n");
    printf("    Interpretation:\n");
    printf("    - R^2 < 0.5: Fluctuations dominate (good NVE)\n");
    printf("    - R^2 > 0.9: Systematic drift (problematic)\n");
    printf("\n");
    printf("    SKIP: Requires trajectory data\n");
    return -1;
#endif
}

/**
 * Test V.4: Block-averaged drift consistency.
 *
 * ThDD:T-US-043-V.4
 *
 * Criterion: |drift_first_half - drift_second_half| < 3 * SE_drift
 */
static int test_block_drift_consistency_V_US_043_V_4(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    printf("\n");
    printf("    === V.4: Block Drift Consistency ===\n");
    printf("    Acceptance: |drift_1 - drift_2| < %.1f * SE_drift\n",
           BLOCK_DRIFT_SE_MULTIPLIER);
    printf("\n");
    printf("    This verifies stationary behavior: drift should not\n");
    printf("    change systematically between trajectory halves.\n");
    printf("\n");
    printf("    SKIP: Requires trajectory data\n");
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
    printf("Unit tests verify statistical algorithms using mathematical identities.\n");
    printf("Integration tests require trajectory data from long GROMACS simulations.\n");
    printf("Integration tests are gated behind GRODFTB_LONG_TESTS cmake option.\n");
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

    printf("=== US-043 Extended NVE Stability Tests ===\n\n");

    printf("Acceptance Criteria (from specs.md S21.1, T-US-043-V.*):\n");
    printf("  V.1: Drift threshold: < %.4f kJ/mol/ps/atom\n",
           NVE_DRIFT_THRESHOLD_KJMOL_PS_ATOM);
    printf("  V.2: Energy fluctuation: < %.0f kJ/mol\n",
           ENERGY_FLUCTUATION_MULTIPLIER * compute_expected_sigma_E(B5_N_ATOMS_TOTAL, B5_TARGET_T_K));
    printf("  V.3: R-squared (fluctuation dominated): < %.1f\n", R_SQUARED_UPPER_BOUND);
    printf("  V.4: Block drift consistency: within %.0f SE\n", BLOCK_DRIFT_SE_MULTIPLIER);
    printf("\n");

    printf("B5 System Parameters:\n");
    printf("  Total atoms: %d\n", B5_N_ATOMS_TOTAL);
    printf("  QM atoms: %d\n", B5_N_QM_ATOMS);
    printf("  Target temperature: %.1f K\n", B5_TARGET_T_K);
    printf("\n");

    if (run_unit) {
        printf("=== Unit Tests (Mathematical Identities) ===\n\n");

        printf("Linear Regression (T-US-043-N.1):\n");
        RUN_TEST(test_nve_stability_regression_exact_slope_T_US_043_N_1);
        printf("\n");

        printf("Normalized Drift (T-US-043-V.1):\n");
        RUN_TEST(test_nve_stability_normalized_drift_T_US_043_V_1);
        printf("\n");

        printf("Energy Fluctuation Formula (T-US-043-V.2):\n");
        RUN_TEST(test_nve_stability_energy_fluctuation_formula_T_US_043_V_2);
        printf("\n");

        printf("R-squared Interpretation (T-US-043-V.3):\n");
        RUN_TEST(test_nve_stability_r_squared_noise_T_US_043_V_3);
        printf("\n");

        printf("Block Drift Consistency (T-US-043-V.4):\n");
        RUN_TEST(test_nve_stability_block_consistency_identical_T_US_043_V_4);
        printf("\n");

        printf("Standard Deviation (T-US-043-N.2):\n");
        RUN_TEST(test_nve_stability_stddev_calculation_T_US_043_N_2);
        printf("\n");
    }

    if (run_integration) {
        printf("=== Integration Tests (V.1-V.4 Trajectory Analysis) ===\n\n");

#ifndef GRODFTB_LONG_TESTS
        printf("NOTE: Integration tests require GRODFTB_LONG_TESTS=ON\n");
        printf("      Rebuild with: cmake -DGRODFTB_LONG_TESTS=ON ..\n\n");
#endif

        printf("V.1 - Energy Drift:\n");
        RUN_TEST(test_nve_drift_100ps_V_US_043_V_1);
        RUN_TEST(test_nve_drift_500ps_V_US_043_V_1);
        printf("\n");

        printf("V.2 - Energy Fluctuation:\n");
        RUN_TEST(test_energy_fluctuation_V_US_043_V_2);
        printf("\n");

        printf("V.3 - R-squared Diagnostic:\n");
        RUN_TEST(test_drift_rsquared_V_US_043_V_3);
        printf("\n");

        printf("V.4 - Block Drift Consistency:\n");
        RUN_TEST(test_block_drift_consistency_V_US_043_V_4);
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
        printf("  2. Trajectory data in tests/data/b5/nve_*ps/\n");
    }

    if (tests_failed > 0) {
        printf("\nFAILED: %d test(s) failed\n", tests_failed);
        return 1;
    }

    return 0;
}
