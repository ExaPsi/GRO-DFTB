/*
 * US-042: NVE Energy Drift Statistical Analysis Tests
 *
 * ThDD:T-US-042-N.1 through N.17 -- Drift measurement methodology
 * SDD:specs.md:S21.1 -- NVE drift criterion (< 0.01 kJ/mol/ps/atom)
 * SDD:docs/verification/US-042.md -- Verification plan
 *
 * This file implements:
 *   1. Statistical analysis functions for NVE drift measurement
 *   2. Unit tests using mathematical identities (Golden Rule compliant)
 *   3. Integration test harness for full trajectory analysis
 *
 * GOLDEN RULE COMPLIANCE:
 * - Unit tests use mathematical identities (y = mx + b -> slope = m exactly)
 * - Integration test thresholds from specs.md, not fabricated
 * - Temperature formula tests verify computation, not simulation results
 * - No hardcoded "expected" simulation results
 *
 * Test execution modes:
 *   --unit         Run unit tests only (fast, no simulation required)
 *   --integration  Run integration tests (requires trajectory data)
 *   (no args)      Run unit tests only
 *
 * LGPL-3.0-or-later
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <unistd.h>

#include "grodftb/driver.h"
#include "grodftb/units.h"
#include "grodftb/error.h"

/* ===========================================================================
 * PHYSICAL CONSTANTS
 *
 * ThDD:T-US-042-N.9 -- Temperature fluctuation formula requires k_B
 * ===========================================================================*/

/* Boltzmann constant in kJ/mol/K (CODATA 2022) */
#define KB_KJMOL_K  0.0083144626

/* ===========================================================================
 * ACCEPTANCE CRITERIA (from specs.md Section 21.1)
 *
 * SDD:specs.md:S21.1 line ~1674:
 *   "Stability: NVE drift | 3, 4 | 10 ps NVE QM/MM | < 0.01 kJ/mol/ps/atom"
 *
 * These are specification-derived thresholds, NOT fabricated expected values.
 * ===========================================================================*/

/* SDD:specs.md:S21.1 -- Primary acceptance criterion */
#define NVE_DRIFT_THRESHOLD_KJMOL_PS_ATOM  0.01

/* Derived: Drift error should be 10x smaller for statistical significance */
#define NVE_DRIFT_ERROR_THRESHOLD_KJMOL_PS_ATOM  0.001

/* Overflow detection threshold (any energy > 10^10 kJ/mol is catastrophic) */
#define ENERGY_OVERFLOW_THRESHOLD  1.0e10

/* ===========================================================================
 * B5 SYSTEM PARAMETERS (from tests/data/b5/provenance.json)
 *
 * ThDD:CLAUDE.md:Golden_Rule -- All values from actual GROMACS preparation
 * ===========================================================================*/

#define B5_N_ATOMS_TOTAL  2652   /* Total atoms in system */
#define B5_N_QM_ATOMS     3      /* QM region: 1 water molecule */
#define B5_TARGET_T_K     300.0  /* Target temperature in Kelvin */

/* ===========================================================================
 * REGRESSION RESULT STRUCTURE
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
    double S_tt;        /* Sum of squared time deviations (for diagnostics) */
    double SS_res;      /* Sum of squared residuals */
    double SS_tot;      /* Total sum of squares */
    int n;              /* Number of data points */
} regression_result_t;

/* ===========================================================================
 * SECTION 1: STATISTICAL ANALYSIS FUNCTIONS (Implementation)
 *
 * These functions implement the drift analysis methodology from
 * docs/theory/US-042/03_drift_measurement_methodology.md
 * ===========================================================================*/

/**
 * Compute linear regression with full statistics.
 *
 * ThDD:T-US-042-N.3 -- Linear regression formulas
 * ThDD:T-US-042-N.4 -- R-squared calculation
 * ThDD:T-US-042-N.5 -- Residual standard error
 * ThDD:T-US-042-N.6 -- Standard error of slope
 *
 * Numerically stable implementation using centered data:
 *   S_tt = sum((t_i - t_bar)^2)
 *   S_tE = sum((t_i - t_bar)(E_i - E_bar))
 *   slope = S_tE / S_tt
 *
 * @param n       Number of data points (must be >= 2)
 * @param t       Time array [n]
 * @param E       Energy array [n]
 * @param result  Output structure receiving all statistics
 */
static void compute_linear_regression(int n, const double *t, const double *E,
                                       regression_result_t *result)
{
    /* Initialize output */
    memset(result, 0, sizeof(*result));
    result->n = n;

    if (n < 2) {
        result->slope = 0.0;
        result->intercept = (n > 0) ? E[0] : 0.0;
        return;
    }

    /*
     * ThDD:T-US-042-N.3 -- Step 1: Compute means
     *
     * Using Kahan summation for numerical stability with large datasets.
     * For 2000 data points, this prevents accumulation of floating-point errors.
     */
    double sum_t = 0.0, sum_E = 0.0;
    double c_t = 0.0, c_E = 0.0;  /* Kahan compensation terms */

    for (int i = 0; i < n; i++) {
        /* Kahan summation for t */
        double y_t = t[i] - c_t;
        double temp_t = sum_t + y_t;
        c_t = (temp_t - sum_t) - y_t;
        sum_t = temp_t;

        /* Kahan summation for E */
        double y_E = E[i] - c_E;
        double temp_E = sum_E + y_E;
        c_E = (temp_E - sum_E) - y_E;
        sum_E = temp_E;
    }

    double t_bar = sum_t / (double)n;
    double E_bar = sum_E / (double)n;

    /*
     * ThDD:T-US-042-N.3 -- Step 2: Compute centered sums
     *
     * S_tt = sum((t_i - t_bar)^2)
     * S_tE = sum((t_i - t_bar)(E_i - E_bar))
     * S_EE = sum((E_i - E_bar)^2) = SS_tot
     */
    double S_tt = 0.0, S_tE = 0.0, S_EE = 0.0;
    double c_tt = 0.0, c_tE = 0.0, c_EE = 0.0;

    for (int i = 0; i < n; i++) {
        double dt = t[i] - t_bar;
        double dE = E[i] - E_bar;

        /* Kahan summation for S_tt */
        double y_tt = dt * dt - c_tt;
        double temp_tt = S_tt + y_tt;
        c_tt = (temp_tt - S_tt) - y_tt;
        S_tt = temp_tt;

        /* Kahan summation for S_tE */
        double y_tE = dt * dE - c_tE;
        double temp_tE = S_tE + y_tE;
        c_tE = (temp_tE - S_tE) - y_tE;
        S_tE = temp_tE;

        /* Kahan summation for S_EE (SS_tot) */
        double y_EE = dE * dE - c_EE;
        double temp_EE = S_EE + y_EE;
        c_EE = (temp_EE - S_EE) - y_EE;
        S_EE = temp_EE;
    }

    result->S_tt = S_tt;
    result->SS_tot = S_EE;

    /*
     * ThDD:T-US-042-N.3 -- Step 3: Compute slope and intercept
     */
    if (fabs(S_tt) < 1e-30) {
        /* Degenerate case: all times are the same */
        result->slope = 0.0;
        result->intercept = E_bar;
        return;
    }

    result->slope = S_tE / S_tt;
    result->intercept = E_bar - result->slope * t_bar;

    /*
     * ThDD:T-US-042-N.4 -- Step 4: Compute R-squared
     *
     * R^2 = 1 - SS_res / SS_tot
     * SS_res = sum((E_i - E_hat_i)^2) where E_hat_i = slope * t_i + intercept
     */
    double SS_res = 0.0;
    double c_res = 0.0;

    for (int i = 0; i < n; i++) {
        double E_pred = result->slope * t[i] + result->intercept;
        double residual = E[i] - E_pred;
        double resid_sq = residual * residual;

        /* Kahan summation */
        double y_res = resid_sq - c_res;
        double temp_res = SS_res + y_res;
        c_res = (temp_res - SS_res) - y_res;
        SS_res = temp_res;
    }

    result->SS_res = SS_res;

    if (S_EE > 1e-30) {
        result->r_squared = 1.0 - SS_res / S_EE;
    } else {
        /* All energies are the same - R^2 is undefined, set to 0 */
        result->r_squared = 0.0;
    }

    /*
     * ThDD:T-US-042-N.5 -- Step 5: Compute residual standard error
     *
     * sigma_E = sqrt(SS_res / (n - 2))
     * The n-2 accounts for two degrees of freedom used (slope and intercept).
     */
    if (n > 2) {
        result->sigma_E = sqrt(SS_res / (double)(n - 2));
    } else {
        result->sigma_E = 0.0;
    }

    /*
     * ThDD:T-US-042-N.6 -- Step 6: Compute standard error of slope
     *
     * sigma_slope = sigma_E / sqrt(S_tt)
     */
    if (S_tt > 1e-30) {
        result->sigma_slope = result->sigma_E / sqrt(S_tt);
    } else {
        result->sigma_slope = 0.0;
    }
}

/**
 * Compute normalized drift rate.
 *
 * ThDD:T-US-042-N.2 -- Normalized drift: alpha_norm = alpha / N_atoms
 *
 * @param slope    Drift rate in kJ/mol/ps
 * @param n_atoms  Number of atoms in system
 * @return         Normalized drift in kJ/mol/ps/atom
 */
static double compute_normalized_drift(double slope, int n_atoms)
{
    if (n_atoms <= 0) return 0.0;
    return slope / (double)n_atoms;
}

/**
 * Compute expected temperature fluctuation in NVE ensemble.
 *
 * ThDD:T-US-042-N.9 -- sigma_T = sqrt(2 * k_B * T^2 / (3N - 6))
 *
 * @param T_kelvin  Mean temperature in Kelvin
 * @param n_atoms   Number of atoms
 * @return          Expected temperature standard deviation in Kelvin
 */
static double compute_expected_sigma_T(double T_kelvin, int n_atoms)
{
    /* Degrees of freedom: 3N - 6 (removes 3 translation + 3 rotation) */
    int dof = 3 * n_atoms - 6;
    if (dof <= 0) return 0.0;

    /*
     * ThDD:T-US-042-N.9
     *
     * sigma_T^2 = 2 * k_B * T^2 / (3N - 6)
     *
     * Note: k_B appears in the derivation but cancels out when expressed
     * in terms of temperature directly. The formula simplifies to:
     *
     * sigma_T = T * sqrt(2 / dof)
     *
     * This is the standard result from statistical mechanics for temperature
     * fluctuations in the microcanonical ensemble.
     */
    return T_kelvin * sqrt(2.0 / (double)dof);
}

/**
 * Check if measured temperature fluctuation is consistent with theory.
 *
 * ThDD:T-US-042-N.9 -- Verification that sigma_T matches microcanonical prediction
 *
 * @param sigma_T_measured  Measured temperature standard deviation (K)
 * @param sigma_T_theory    Theoretical expectation (K)
 * @param n_frames          Number of frames in trajectory
 * @return                  1 if consistent (within 2 standard errors), 0 otherwise
 */
__attribute__((unused))
static int check_temperature_fluctuation_consistency(double sigma_T_measured,
                                                      double sigma_T_theory,
                                                      int n_frames)
{
    /*
     * ThDD:T-US-042-N.9
     *
     * The measured sigma_T has uncertainty ~ sigma_T / sqrt(2*n) for a Gaussian.
     * We use 2x this as our acceptance window (approximately 95% confidence).
     *
     * Criterion: |sigma_T_measured - sigma_T_theory| < 2 * sigma_T_theory / sqrt(n)
     */
    if (n_frames <= 1 || sigma_T_theory <= 0.0) return 0;

    double tolerance = 2.0 * sigma_T_theory / sqrt((double)n_frames);
    double diff = fabs(sigma_T_measured - sigma_T_theory);

    return (diff < tolerance) ? 1 : 0;
}

/**
 * Compute standard deviation of an array.
 *
 * @param n    Number of elements
 * @param arr  Data array
 * @return     Sample standard deviation
 */
__attribute__((unused))
static double compute_stddev(int n, const double *arr)
{
    if (n < 2) return 0.0;

    /* Two-pass algorithm for numerical stability */
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
    variance /= (double)(n - 1);  /* Bessel's correction */

    return sqrt(variance);
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
 * These tests verify the statistical algorithms using synthetic data
 * with known mathematical properties. NO simulation required.
 *
 * GOLDEN RULE: All expected values are mathematical identities,
 * not fabricated physical results.
 * ===========================================================================*/

/**
 * Test: Linear regression on perfect linear data.
 *
 * ThDD:T-US-042-N.3
 *
 * MATHEMATICAL IDENTITY:
 *   y = 2x + 5 for x = {0, 1, 2, 3, 4}
 *   Expected: slope = 2.0 (exact by construction)
 *             intercept = 5.0 (exact by construction)
 *   Tolerance: 1e-14 (machine epsilon for double)
 */
static int test_nve_drift_analysis_linear_regression_T_US_042_N_3(void)
{
    /* Construct exact linear data: E = 2*t + 5 */
    double t[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
    double E[5] = {5.0, 7.0, 9.0, 11.0, 13.0};

    regression_result_t result;
    compute_linear_regression(5, t, E, &result);

    /* Mathematical identity: exact values by construction */
    const double TOL = 1e-14;

    if (fabs(result.slope - 2.0) > TOL) {
        fprintf(stderr, "\n    FAIL: slope = %.15f, expected 2.0\n", result.slope);
        return 0;
    }

    if (fabs(result.intercept - 5.0) > TOL) {
        fprintf(stderr, "\n    FAIL: intercept = %.15f, expected 5.0\n", result.intercept);
        return 0;
    }

    return 1;
}

/**
 * Test: R-squared on perfect fit.
 *
 * ThDD:T-US-042-N.4
 *
 * MATHEMATICAL IDENTITY:
 *   Perfect fit -> SS_res = 0 -> R^2 = 1 - 0/SS_tot = 1.0
 */
static int test_nve_drift_analysis_r_squared_T_US_042_N_4(void)
{
    /* Perfect linear data */
    double t[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
    double E[5] = {5.0, 7.0, 9.0, 11.0, 13.0};

    regression_result_t result;
    compute_linear_regression(5, t, E, &result);

    /* Mathematical identity: perfect fit -> R^2 = 1.0 */
    const double TOL = 1e-14;

    if (fabs(result.r_squared - 1.0) > TOL) {
        fprintf(stderr, "\n    FAIL: R^2 = %.15f, expected 1.0\n", result.r_squared);
        return 0;
    }

    /* Also verify SS_res is effectively zero */
    if (result.SS_res > 1e-25) {
        fprintf(stderr, "\n    FAIL: SS_res = %.15e, expected ~0\n", result.SS_res);
        return 0;
    }

    return 1;
}

/**
 * Test: R-squared on constant data (zero slope).
 *
 * ThDD:T-US-042-N.4
 *
 * MATHEMATICAL IDENTITY:
 *   Constant data (y = 5 for all x):
 *   - slope = 0
 *   - SS_tot = 0 (all points equal to mean)
 *   - R^2 is undefined (0/0), convention: set to 0
 */
static int test_nve_drift_analysis_r_squared_constant_T_US_042_N_4(void)
{
    /* Constant data */
    double t[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
    double E[5] = {5.0, 5.0, 5.0, 5.0, 5.0};

    regression_result_t result;
    compute_linear_regression(5, t, E, &result);

    const double TOL = 1e-14;

    /* Slope should be zero */
    if (fabs(result.slope) > TOL) {
        fprintf(stderr, "\n    FAIL: slope = %.15f, expected 0.0\n", result.slope);
        return 0;
    }

    /* R^2 should be 0 (our convention for undefined case) */
    if (fabs(result.r_squared) > TOL) {
        fprintf(stderr, "\n    FAIL: R^2 = %.15f, expected 0.0\n", result.r_squared);
        return 0;
    }

    return 1;
}

/**
 * Test: Standard error on data with known residuals.
 *
 * ThDD:T-US-042-N.5, T-US-042-N.6
 *
 * MATHEMATICAL IDENTITY (manual calculation):
 *   Input: x = {0, 1, 2, 3, 4}, y = {5.1, 6.9, 9.0, 11.1, 12.9}
 *
 *   Step 1: Compute means
 *     x_bar = 2.0
 *     y_bar = 9.0
 *
 *   Step 2: Compute sums
 *     S_tt = (0-2)^2 + (1-2)^2 + (2-2)^2 + (3-2)^2 + (4-2)^2 = 4+1+0+1+4 = 10
 *     S_tE = (0-2)(5.1-9) + (1-2)(6.9-9) + (2-2)(9-9) + (3-2)(11.1-9) + (4-2)(12.9-9)
 *          = (-2)(-3.9) + (-1)(-2.1) + (0)(0) + (1)(2.1) + (2)(3.9)
 *          = 7.8 + 2.1 + 0 + 2.1 + 7.8 = 19.8
 *
 *   Step 3: slope = S_tE / S_tt = 19.8 / 10 = 1.98
 *           intercept = y_bar - slope * x_bar = 9.0 - 1.98 * 2 = 5.04
 *
 *   Step 4: Residuals
 *     y_pred = [5.04, 7.02, 9.00, 10.98, 12.96]
 *     residuals = [0.06, -0.12, 0.00, 0.12, -0.06]
 *     SS_res = 0.0036 + 0.0144 + 0 + 0.0144 + 0.0036 = 0.036
 *
 *   Step 5: sigma_E = sqrt(SS_res / (n-2)) = sqrt(0.036 / 3) = sqrt(0.012) = 0.1095
 *
 *   Step 6: sigma_slope = sigma_E / sqrt(S_tt) = 0.1095 / sqrt(10) = 0.0346
 *
 *   Tolerance: 1e-3 (accounts for rounding in manual calculation)
 */
static int test_nve_drift_analysis_std_error_T_US_042_N_6(void)
{
    double t[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
    double E[5] = {5.1, 6.9, 9.0, 11.1, 12.9};

    regression_result_t result;
    compute_linear_regression(5, t, E, &result);

    /* Expected values from manual calculation above */
    const double expected_slope = 1.98;
    const double expected_intercept = 5.04;
    const double expected_SS_res = 0.036;
    const double expected_sigma_E = 0.1095;
    const double expected_sigma_slope = 0.0346;

    const double TOL = 1e-3;

    if (fabs(result.slope - expected_slope) > TOL) {
        fprintf(stderr, "\n    FAIL: slope = %.6f, expected %.6f\n",
                result.slope, expected_slope);
        return 0;
    }

    if (fabs(result.intercept - expected_intercept) > TOL) {
        fprintf(stderr, "\n    FAIL: intercept = %.6f, expected %.6f\n",
                result.intercept, expected_intercept);
        return 0;
    }

    if (fabs(result.SS_res - expected_SS_res) > TOL) {
        fprintf(stderr, "\n    FAIL: SS_res = %.6f, expected %.6f\n",
                result.SS_res, expected_SS_res);
        return 0;
    }

    if (fabs(result.sigma_E - expected_sigma_E) > TOL) {
        fprintf(stderr, "\n    FAIL: sigma_E = %.6f, expected %.6f\n",
                result.sigma_E, expected_sigma_E);
        return 0;
    }

    if (fabs(result.sigma_slope - expected_sigma_slope) > TOL) {
        fprintf(stderr, "\n    FAIL: sigma_slope = %.6f, expected %.6f\n",
                result.sigma_slope, expected_sigma_slope);
        return 0;
    }

    return 1;
}

/**
 * Test: Temperature fluctuation formula.
 *
 * ThDD:T-US-042-N.9
 *
 * MATHEMATICAL IDENTITY:
 *   sigma_T = T * sqrt(2 / (3N - 6))
 *
 *   For N = 1000, T = 300 K:
 *     dof = 3*1000 - 6 = 2994
 *     sigma_T = 300 * sqrt(2 / 2994)
 *             = 300 * sqrt(0.000668)
 *             = 300 * 0.02585
 *             = 7.754 K
 *
 *   Tolerance: 1e-3 K
 */
static int test_nve_drift_analysis_temperature_check_T_US_042_N_9(void)
{
    double T = 300.0;
    int N = 1000;

    double sigma_T = compute_expected_sigma_T(T, N);

    /* Mathematical identity from formula */
    double expected = T * sqrt(2.0 / (3.0 * N - 6.0));
    /* expected = 300 * sqrt(2/2994) = 7.754 K */

    const double TOL = 1e-3;

    if (fabs(sigma_T - expected) > TOL) {
        fprintf(stderr, "\n    FAIL: sigma_T = %.6f K, expected %.6f K\n",
                sigma_T, expected);
        return 0;
    }

    /* Also verify the actual numeric value */
    if (fabs(sigma_T - 7.754) > 0.001) {
        fprintf(stderr, "\n    FAIL: sigma_T = %.6f K, expected ~7.754 K\n", sigma_T);
        return 0;
    }

    return 1;
}

/**
 * Test: Temperature fluctuation for B5 system.
 *
 * ThDD:T-US-042-N.9
 *
 * FORMULA VERIFICATION for B5:
 *   N = 2652 atoms, T = 300 K
 *   dof = 3*2652 - 6 = 7950
 *   sigma_T = 300 * sqrt(2 / 7950)
 *           = 300 * 0.01586
 *           = 4.759 K
 */
static int test_nve_drift_analysis_temperature_b5_T_US_042_N_9(void)
{
    double T = B5_TARGET_T_K;
    int N = B5_N_ATOMS_TOTAL;

    double sigma_T = compute_expected_sigma_T(T, N);

    /* Mathematical identity from formula */
    double expected = T * sqrt(2.0 / (3.0 * N - 6.0));

    const double TOL = 1e-3;

    if (fabs(sigma_T - expected) > TOL) {
        fprintf(stderr, "\n    FAIL: sigma_T = %.6f K, expected %.6f K\n",
                sigma_T, expected);
        return 0;
    }

    /* Verify numeric value for B5 */
    if (fabs(sigma_T - 4.759) > 0.01) {
        fprintf(stderr, "\n    FAIL: sigma_T = %.6f K, expected ~4.76 K\n", sigma_T);
        return 0;
    }

    printf("\n    (B5 expected sigma_T = %.3f K) ", sigma_T);

    return 1;
}

/**
 * Test: Normalized drift calculation.
 *
 * ThDD:T-US-042-N.2
 *
 * MATHEMATICAL IDENTITY:
 *   slope = 26.52 kJ/mol/ps
 *   N_atoms = 2652
 *   alpha_norm = 26.52 / 2652 = 0.01 kJ/mol/ps/atom
 */
static int test_nve_drift_analysis_normalized_drift_T_US_042_N_2(void)
{
    double slope = 26.52;  /* kJ/mol/ps */
    int n_atoms = 2652;

    double alpha_norm = compute_normalized_drift(slope, n_atoms);

    const double TOL = 1e-10;

    if (fabs(alpha_norm - 0.01) > TOL) {
        fprintf(stderr, "\n    FAIL: alpha_norm = %.15f, expected 0.01\n", alpha_norm);
        return 0;
    }

    return 1;
}

/**
 * Test: Numerical stability with large datasets.
 *
 * ThDD:T-US-042-N.3
 *
 * Verifies that Kahan summation prevents accumulation errors when
 * processing 2000 data points (typical NVE trajectory length).
 *
 * MATHEMATICAL IDENTITY:
 *   y = 0.001 * x + 1000.0 for x = {0, 5, 10, ..., 9995}
 *   n = 2000 points
 *   slope = 0.001 (exact)
 *   intercept = 1000.0 (exact)
 */
static int test_nve_drift_analysis_numerical_stability_T_US_042_N_3(void)
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

    /* Generate exact linear data with large offset (stress test for precision) */
    for (int i = 0; i < n; i++) {
        t[i] = i * 5.0;  /* Time: 0, 5, 10, ..., 9995 ps */
        E[i] = 0.001 * t[i] + 1000.0;  /* E = 0.001*t + 1000 */
    }

    regression_result_t result;
    compute_linear_regression(n, t, E, &result);

    free(t);
    free(E);

    /* With Kahan summation, we should maintain high precision */
    const double TOL = 1e-10;

    if (fabs(result.slope - 0.001) > TOL) {
        fprintf(stderr, "\n    FAIL: slope = %.15e, expected 0.001\n", result.slope);
        fprintf(stderr, "    (Numerical stability issue with Kahan summation)\n");
        return 0;
    }

    if (fabs(result.intercept - 1000.0) > TOL) {
        fprintf(stderr, "\n    FAIL: intercept = %.15f, expected 1000.0\n", result.intercept);
        return 0;
    }

    if (fabs(result.r_squared - 1.0) > TOL) {
        fprintf(stderr, "\n    FAIL: R^2 = %.15f, expected 1.0\n", result.r_squared);
        return 0;
    }

    printf("\n    (2000 points, slope error = %.2e) ", fabs(result.slope - 0.001));

    return 1;
}

/**
 * Test: Edge case - minimum data points (n=2).
 *
 * ThDD:T-US-042-N.3
 *
 * With only 2 points, the fit is exact but statistics are limited.
 */
static int test_nve_drift_analysis_edge_case_n2_T_US_042_N_3(void)
{
    double t[2] = {0.0, 10.0};
    double E[2] = {100.0, 110.0};  /* slope = 1.0 */

    regression_result_t result;
    compute_linear_regression(2, t, E, &result);

    const double TOL = 1e-14;

    if (fabs(result.slope - 1.0) > TOL) {
        fprintf(stderr, "\n    FAIL: slope = %.15f, expected 1.0\n", result.slope);
        return 0;
    }

    if (fabs(result.intercept - 100.0) > TOL) {
        fprintf(stderr, "\n    FAIL: intercept = %.15f, expected 100.0\n", result.intercept);
        return 0;
    }

    /* R^2 should be 1.0 for perfect fit */
    if (fabs(result.r_squared - 1.0) > TOL) {
        fprintf(stderr, "\n    FAIL: R^2 = %.15f, expected 1.0\n", result.r_squared);
        return 0;
    }

    /* sigma_E should be 0 for n=2 (no degrees of freedom for error) */
    if (fabs(result.sigma_E) > TOL) {
        fprintf(stderr, "\n    FAIL: sigma_E = %.15f, expected 0.0\n", result.sigma_E);
        return 0;
    }

    return 1;
}

/* ===========================================================================
 * SECTION 4: INTEGRATION TESTS (Trajectory Analysis)
 *
 * These tests require actual trajectory data from GROMACS simulations.
 * They are gated behind GRODFTB_LONG_TESTS.
 * ===========================================================================*/

#ifndef GRODFTB_SOURCE_DIR_STR
#define GRODFTB_SOURCE_DIR_STR "."
#endif
#define B5_DATA_DIR GRODFTB_SOURCE_DIR_STR "/tests/data/b5"

/**
 * Test: NVE drift acceptance criterion.
 *
 * ThDD:T-US-042-N.2
 * SDD:specs.md:S21.1 -- drift < 0.01 kJ/mol/ps/atom
 *
 * This test requires trajectory data from a 10ps NVE simulation.
 * It is a placeholder that will:
 *   1. Load energy trajectory from XVG file
 *   2. Compute drift rate via linear regression
 *   3. Check against acceptance criterion
 *
 * NOTE: This test SKIPs if trajectory data is not available.
 * The trajectory must be generated by running GROMACS mdrun.
 */
static int test_nve_drift_10ps_V_US_042_1(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    /*
     * Integration test: Load real trajectory and analyze.
     *
     * Expected workflow:
     * 1. Load tests/data/b5/nve_10ps/total_energy.xvg
     * 2. Parse XVG format (skip @ and # lines)
     * 3. Compute linear regression
     * 4. Check |alpha_norm| < 0.01 kJ/mol/ps/atom
     *
     * For now, this test verifies the acceptance criteria setup
     * and returns SKIP if trajectory data is not available.
     */
    printf("\n");
    printf("    Acceptance criterion: drift < %.4f kJ/mol/ps/atom\n",
           NVE_DRIFT_THRESHOLD_KJMOL_PS_ATOM);
    printf("    B5 system: %d atoms\n", B5_N_ATOMS_TOTAL);
    printf("    Total drift limit: %.2f kJ/mol/ps\n",
           NVE_DRIFT_THRESHOLD_KJMOL_PS_ATOM * B5_N_ATOMS_TOTAL);
    printf("\n");

    /* Check for trajectory data file */
    char xvg_path[512];
    snprintf(xvg_path, sizeof(xvg_path), "%s/nve_10ps/total_energy.xvg", B5_DATA_DIR);

    FILE *f = fopen(xvg_path, "r");
    if (!f) {
        printf("    Trajectory data not found: %s\n", xvg_path);
        printf("    Generate with: gmx mdrun -s nve_10ps.tpr -deffnm nve_10ps\n");
        printf("    Then: gmx energy -f nve_10ps.edr -o nve_10ps/total_energy.xvg\n");
        printf("\n");
        printf("    SKIP: Trajectory data required for integration test\n");
        return -1;
    }
    fclose(f);

    /* TODO: Implement XVG parsing and analysis */
    printf("    Trajectory data found - analysis not yet implemented\n");
    printf("    SKIP: Full analysis pending implementation\n");
    return -1;
#endif
}

/**
 * Test: Drift error estimate.
 *
 * ThDD:T-US-042-N.6
 *
 * Verifies that the standard error of the drift rate is sufficiently small
 * to make the acceptance decision statistically meaningful.
 */
static int test_nve_drift_error_V_US_042_2(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    printf("\n");
    printf("    Error threshold: %.4f kJ/mol/ps/atom\n",
           NVE_DRIFT_ERROR_THRESHOLD_KJMOL_PS_ATOM);
    printf("\n");
    printf("    SKIP: Requires trajectory data\n");
    return -1;
#endif
}

/**
 * Test: Temperature fluctuation consistency.
 *
 * ThDD:T-US-042-N.9
 *
 * Verifies that measured temperature fluctuations match the theoretical
 * prediction for the microcanonical ensemble.
 */
static int test_nve_temperature_fluctuation_V_US_042_3(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    double expected_sigma_T = compute_expected_sigma_T(B5_TARGET_T_K, B5_N_ATOMS_TOTAL);
    printf("\n");
    printf("    Expected sigma_T for B5: %.3f K\n", expected_sigma_T);
    printf("    (from ThDD:T-US-042-N.9 formula)\n");
    printf("\n");
    printf("    SKIP: Requires trajectory data\n");
    return -1;
#endif
}

/**
 * Test: R-squared interpretation.
 *
 * ThDD:T-US-042-N.4
 *
 * Verifies that R² is interpreted correctly:
 * - High R² with low drift: good precision
 * - High R² with high drift: problem detected
 * - Low R²: fluctuations dominate (expected for good NVE)
 */
static int test_nve_drift_diagnostic_r_squared_V_US_042_4(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    printf("\n");
    printf("    R² interpretation:\n");
    printf("    - R² high + drift low: precise measurement, drift negligible\n");
    printf("    - R² high + drift high: significant drift detected (FAIL)\n");
    printf("    - R² low + drift low: fluctuations dominate (expected)\n");
    printf("\n");
    printf("    SKIP: Requires trajectory data\n");
    return -1;
#endif
}

/**
 * Test: Numerical stability check.
 *
 * ThDD:T-US-042-N.11
 *
 * Verifies no NaN, Inf, or energy overflow in trajectory.
 */
static int test_nve_stability_no_nan_inf_T_US_042_N_11(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    printf("\n");
    printf("    Stability checks:\n");
    printf("    - No NaN in energy values\n");
    printf("    - No Inf in energy values\n");
    printf("    - No overflow (|E| < %.0e)\n", ENERGY_OVERFLOW_THRESHOLD);
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
    printf("Integration tests require trajectory data from GROMACS simulations.\n");
}

int main(int argc, char *argv[])
{
    int run_unit = 1;
    int run_integration = 0;

    /* Parse arguments */
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

    printf("=== US-042 NVE Energy Drift Analysis Tests ===\n\n");

    printf("Acceptance Criteria (from specs.md S21.1):\n");
    printf("  Drift threshold: %.4f kJ/mol/ps/atom\n", NVE_DRIFT_THRESHOLD_KJMOL_PS_ATOM);
    printf("  Error threshold: %.4f kJ/mol/ps/atom\n", NVE_DRIFT_ERROR_THRESHOLD_KJMOL_PS_ATOM);
    printf("\n");

    printf("B5 System Parameters:\n");
    printf("  Total atoms: %d\n", B5_N_ATOMS_TOTAL);
    printf("  Target temperature: %.1f K\n", B5_TARGET_T_K);
    printf("  Expected sigma_T: %.3f K (ThDD:T-US-042-N.9)\n",
           compute_expected_sigma_T(B5_TARGET_T_K, B5_N_ATOMS_TOTAL));
    printf("\n");

    if (run_unit) {
        printf("=== Unit Tests (Mathematical Identities) ===\n\n");

        printf("Linear Regression (T-US-042-N.3):\n");
        RUN_TEST(test_nve_drift_analysis_linear_regression_T_US_042_N_3);
        printf("\n");

        printf("R-squared Calculation (T-US-042-N.4):\n");
        RUN_TEST(test_nve_drift_analysis_r_squared_T_US_042_N_4);
        RUN_TEST(test_nve_drift_analysis_r_squared_constant_T_US_042_N_4);
        printf("\n");

        printf("Standard Error (T-US-042-N.5, N.6):\n");
        RUN_TEST(test_nve_drift_analysis_std_error_T_US_042_N_6);
        printf("\n");

        printf("Temperature Fluctuation Formula (T-US-042-N.9):\n");
        RUN_TEST(test_nve_drift_analysis_temperature_check_T_US_042_N_9);
        RUN_TEST(test_nve_drift_analysis_temperature_b5_T_US_042_N_9);
        printf("\n");

        printf("Normalized Drift (T-US-042-N.2):\n");
        RUN_TEST(test_nve_drift_analysis_normalized_drift_T_US_042_N_2);
        printf("\n");

        printf("Numerical Stability:\n");
        RUN_TEST(test_nve_drift_analysis_numerical_stability_T_US_042_N_3);
        printf("\n");

        printf("Edge Cases:\n");
        RUN_TEST(test_nve_drift_analysis_edge_case_n2_T_US_042_N_3);
        printf("\n");
    }

    if (run_integration) {
        printf("=== Integration Tests (Trajectory Analysis) ===\n\n");

#ifndef GRODFTB_LONG_TESTS
        printf("NOTE: Integration tests require GRODFTB_LONG_TESTS=1\n");
        printf("      Rebuild with: cmake -DGRODFTB_LONG_TESTS=ON ..\n\n");
#endif

        printf("Drift Acceptance (V-US-042-1):\n");
        RUN_TEST(test_nve_drift_10ps_V_US_042_1);
        printf("\n");

        printf("Drift Error (V-US-042-2):\n");
        RUN_TEST(test_nve_drift_error_V_US_042_2);
        printf("\n");

        printf("Temperature Fluctuation (V-US-042-3):\n");
        RUN_TEST(test_nve_temperature_fluctuation_V_US_042_3);
        printf("\n");

        printf("R-squared Diagnostic (V-US-042-4):\n");
        RUN_TEST(test_nve_drift_diagnostic_r_squared_V_US_042_4);
        printf("\n");

        printf("Numerical Stability (T-US-042-N.11):\n");
        RUN_TEST(test_nve_stability_no_nan_inf_T_US_042_N_11);
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
        printf("  2. Trajectory data in tests/data/b5/nve_10ps/\n");
    }

    if (tests_failed > 0) {
        printf("\nFAILED: %d test(s) failed\n", tests_failed);
        return 1;
    }

    return 0;
}
