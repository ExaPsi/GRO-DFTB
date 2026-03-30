/*
 * US-043: Temperature Stability Tests (V.8-V.9)
 *
 * ThDD:T-US-043-V.8 -- Temperature mean criterion
 * ThDD:T-US-043-V.9 -- Temperature fluctuation criterion
 * SDD:specs.md:S21.1 -- NVE simulation requirements
 * SDD:docs/verification/US-043.md -- Verification plan
 *
 * This file implements:
 *   1. V.8: Temperature mean (should equal equilibration target)
 *   2. V.9: Temperature fluctuation (consistent with NVE ensemble)
 *
 * GOLDEN RULE COMPLIANCE:
 * - Unit tests use mathematical identities (exact formulas)
 * - Temperature target (300 K) from simulation setup, not fabricated
 * - Fluctuation bounds from statistical mechanics theory
 * - No hardcoded "expected" simulation results
 *
 * LGPL-3.0-or-later
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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
 * ACCEPTANCE CRITERIA (from T-US-043-V.8 and V.9)
 *
 * These are specification-derived thresholds, NOT fabricated expected values.
 * ===========================================================================*/

/* V.8: Temperature mean tolerance (+/- 5 K from target) */
#define TEMPERATURE_MEAN_TOL_K  5.0

/* V.9: Temperature fluctuation upper bound (K) */
#define TEMPERATURE_FLUCTUATION_MAX_K  10.0

/* Target temperature from NVT equilibration (simulation input parameter) */
#define TARGET_TEMPERATURE_K  300.0

/* ===========================================================================
 * B5 SYSTEM PARAMETERS
 * ===========================================================================*/

#define B5_N_ATOMS_TOTAL  2652
#define B5_N_QM_ATOMS     3

/* ===========================================================================
 * HELPER FUNCTIONS
 * ===========================================================================*/

/**
 * Compute expected temperature fluctuation in NVE ensemble.
 *
 * ThDD:T-US-043-V.9 / T-US-042-N.9 -- Microcanonical fluctuation formula
 *
 *   sigma_T = T * sqrt(2 / N_dof)
 *
 * where N_dof = 3*N - 6 for a non-periodic system (3 translation + 3 rotation removed)
 * or N_dof = 3*N - 3 for periodic (only 3 translation removed).
 *
 * @param T_kelvin  Mean temperature in Kelvin
 * @param n_dof     Number of degrees of freedom
 * @return          Expected temperature standard deviation in Kelvin
 */
static double compute_expected_sigma_T(double T_kelvin, int n_dof)
{
    if (n_dof <= 0) return 0.0;

    /*
     * ThDD:T-US-043-V.9
     *
     * From statistical mechanics of the microcanonical ensemble:
     *   <(delta T)^2> = 2 * T^2 / N_dof
     *   sigma_T = T * sqrt(2 / N_dof)
     */
    return T_kelvin * sqrt(2.0 / (double)n_dof);
}

/**
 * Compute degrees of freedom for a periodic system.
 *
 * ThDD:T-US-043-V.9 -- DOF calculation
 *
 * For periodic boundary conditions: N_dof = 3*N - 3
 * (removes 3 translational DOF, rotation is not constrained in PBC)
 *
 * @param n_atoms  Number of atoms
 * @return         Degrees of freedom
 */
static int compute_dof_periodic(int n_atoms)
{
    return 3 * n_atoms - 3;
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
static double compute_stddev(int n, const double *arr)
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

/* ===========================================================================
 * TEST FRAMEWORK
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
 * SECTION 1: UNIT TESTS (Mathematical Identities)
 *
 * GOLDEN RULE: All expected values are mathematical identities,
 * not fabricated physical results.
 * ===========================================================================*/

/**
 * Test: Mean calculation.
 *
 * MATHEMATICAL IDENTITY:
 *   mean({1, 2, 3, 4, 5}) = 15/5 = 3.0
 */
static int test_temperature_mean_calculation_T_US_043_N_1(void)
{
    double data[5] = {1.0, 2.0, 3.0, 4.0, 5.0};

    double mean = compute_mean(5, data);
    double expected = 3.0;

    const double TOL = 1e-14;

    if (fabs(mean - expected) > TOL) {
        fprintf(stderr, "\n    FAIL: mean = %.15f, expected %.15f\n", mean, expected);
        return 0;
    }

    return 1;
}

/**
 * Test: Standard deviation calculation.
 *
 * MATHEMATICAL IDENTITY:
 *   For data {1, 2, 3, 4, 5}:
 *   variance = ((1-3)^2 + (2-3)^2 + (3-3)^2 + (4-3)^2 + (5-3)^2) / 4
 *            = (4 + 1 + 0 + 1 + 4) / 4 = 2.5
 *   stddev = sqrt(2.5) = 1.5811388...
 */
static int test_temperature_stddev_calculation_T_US_043_N_2(void)
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

/**
 * Test: Degrees of freedom for periodic system.
 *
 * ThDD:T-US-043-V.9
 *
 * MATHEMATICAL IDENTITY:
 *   N_dof = 3*N - 3 for periodic system
 *   For B5 (N=2652): N_dof = 3*2652 - 3 = 7953
 */
static int test_temperature_dof_periodic_T_US_043_V_9(void)
{
    int dof = compute_dof_periodic(B5_N_ATOMS_TOTAL);
    int expected = 3 * B5_N_ATOMS_TOTAL - 3;  /* 7953 */

    if (dof != expected) {
        fprintf(stderr, "\n    FAIL: dof = %d, expected %d\n", dof, expected);
        return 0;
    }

    return 1;
}

/**
 * Test: Temperature fluctuation formula.
 *
 * ThDD:T-US-043-V.9
 *
 * MATHEMATICAL IDENTITY:
 *   sigma_T = T * sqrt(2 / N_dof)
 *
 *   For B5 (N_dof = 7953), T = 300 K:
 *     sigma_T = 300 * sqrt(2 / 7953)
 *             = 300 * sqrt(0.0002515)
 *             = 300 * 0.01586
 *             = 4.758 K
 */
static int test_temperature_fluctuation_formula_T_US_043_V_9(void)
{
    int n_dof = compute_dof_periodic(B5_N_ATOMS_TOTAL);
    double sigma_T = compute_expected_sigma_T(TARGET_TEMPERATURE_K, n_dof);

    /* Expected from formula */
    double expected = TARGET_TEMPERATURE_K * sqrt(2.0 / (double)n_dof);

    const double TOL = 1e-10;

    if (fabs(sigma_T - expected) > TOL) {
        fprintf(stderr, "\n    FAIL: sigma_T = %.15f K, expected %.15f K\n",
                sigma_T, expected);
        return 0;
    }

    /* Verify numeric value (~4.76 K) */
    if (sigma_T < 4.5 || sigma_T > 5.0) {
        fprintf(stderr, "\n    FAIL: sigma_T = %.3f K, expected ~4.76 K\n", sigma_T);
        return 0;
    }

    printf("\n    (B5 expected sigma_T = %.3f K) ", sigma_T);
    return 1;
}

/**
 * Test: Temperature fluctuation scaling with system size.
 *
 * ThDD:T-US-043-V.9
 *
 * MATHEMATICAL IDENTITY:
 *   sigma_T proportional to 1/sqrt(N)
 *   If N1 = 4*N2, then sigma_T(N1) = sigma_T(N2) / 2
 */
static int test_temperature_fluctuation_scaling_T_US_043_V_9(void)
{
    double T = 300.0;
    int n1 = 1000;
    int n2 = 4000;  /* 4x larger */

    int dof1 = 3 * n1 - 3;
    int dof2 = 3 * n2 - 3;

    double sigma_T1 = compute_expected_sigma_T(T, dof1);
    double sigma_T2 = compute_expected_sigma_T(T, dof2);

    /* sigma_T2 should be approximately half of sigma_T1 */
    /* Exact ratio: sqrt(dof1/dof2) = sqrt(2997/11997) = sqrt(0.2498) = 0.4998 */
    double expected_ratio = sqrt((double)dof1 / (double)dof2);
    double actual_ratio = sigma_T2 / sigma_T1;

    const double TOL = 1e-6;

    if (fabs(actual_ratio - expected_ratio) > TOL) {
        fprintf(stderr, "\n    FAIL: ratio = %.6f, expected %.6f\n",
                actual_ratio, expected_ratio);
        return 0;
    }

    printf("\n    (scaling: sigma_T(4000)/sigma_T(1000) = %.4f ~ 0.5) ", actual_ratio);
    return 1;
}

/**
 * Test: Temperature bounds are reasonable.
 *
 * ThDD:T-US-043-V.8, V.9
 *
 * Verifies that acceptance thresholds are physically sensible:
 * - Mean tolerance (5 K) >> expected fluctuation (~5 K)
 * - Max fluctuation (10 K) ~ 2 * expected fluctuation
 */
static int test_temperature_bounds_sensible_T_US_043_V_8(void)
{
    int n_dof = compute_dof_periodic(B5_N_ATOMS_TOTAL);
    double expected_sigma_T = compute_expected_sigma_T(TARGET_TEMPERATURE_K, n_dof);

    /* Mean tolerance should be comparable to or larger than fluctuation */
    if (TEMPERATURE_MEAN_TOL_K < expected_sigma_T) {
        fprintf(stderr, "\n    FAIL: Mean tolerance (%.1f K) < expected sigma_T (%.1f K)\n",
                TEMPERATURE_MEAN_TOL_K, expected_sigma_T);
        fprintf(stderr, "    This would cause false failures due to normal fluctuations\n");
        return 0;
    }

    /* Max fluctuation should be about 2x expected (allow for sampling error) */
    if (TEMPERATURE_FLUCTUATION_MAX_K < 1.5 * expected_sigma_T) {
        fprintf(stderr, "\n    FAIL: Max fluctuation (%.1f K) < 1.5 * expected (%.1f K)\n",
                TEMPERATURE_FLUCTUATION_MAX_K, 1.5 * expected_sigma_T);
        fprintf(stderr, "    Threshold too strict for statistical variation\n");
        return 0;
    }

    printf("\n    (tol=%.0f K, max_fluct=%.0f K, expected_sigma=%.1f K) ",
           TEMPERATURE_MEAN_TOL_K, TEMPERATURE_FLUCTUATION_MAX_K, expected_sigma_T);
    return 1;
}

/* ===========================================================================
 * SECTION 2: INTEGRATION TESTS (Trajectory Analysis)
 *
 * These tests require actual trajectory data from GROMACS simulations.
 * Gated behind GRODFTB_LONG_TESTS.
 * ===========================================================================*/

#ifndef GRODFTB_SOURCE_DIR_STR
#define GRODFTB_SOURCE_DIR_STR "."
#endif
#define B5_DATA_DIR GRODFTB_SOURCE_DIR_STR "/tests/data/b5"

/**
 * Test V.8: Temperature mean.
 *
 * ThDD:T-US-043-V.8
 *
 * Criterion: |<T> - T_target| < 5 K
 *
 * NOTE: T_target = 300 K comes from the NVT equilibration setup,
 * not a fabricated expected value.
 */
static int test_temperature_mean_V_US_043_V_8(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    printf("\n");
    printf("    === V.8: Temperature Mean ===\n");
    printf("    Target: %.1f K (from NVT equilibration)\n", TARGET_TEMPERATURE_K);
    printf("    Acceptance: |<T> - %.1f| < %.1f K\n",
           TARGET_TEMPERATURE_K, TEMPERATURE_MEAN_TOL_K);
    printf("\n");

    char xvg_path[512];
    snprintf(xvg_path, sizeof(xvg_path), "%s/nve_500ps/temperature.xvg", B5_DATA_DIR);

    FILE *f = fopen(xvg_path, "r");
    if (!f) {
        printf("    Temperature data not found: %s\n", xvg_path);
        printf("    Generate with: gmx energy -f nve_500ps.edr -o temperature.xvg\n");
        printf("    SKIP: Temperature data required\n");
        return -1;
    }
    fclose(f);

    printf("    Temperature data found - analysis not yet implemented\n");
    printf("    SKIP: Full analysis pending implementation\n");
    return -1;
#endif
}

/**
 * Test V.9: Temperature fluctuation.
 *
 * ThDD:T-US-043-V.9
 *
 * Criterion: sigma_T < 10 K
 *
 * NOTE: Expected fluctuation for B5 at 300 K is ~4.8 K (from formula).
 * The 10 K bound allows for 2x margin.
 */
static int test_temperature_fluctuation_V_US_043_V_9(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    int n_dof = compute_dof_periodic(B5_N_ATOMS_TOTAL);
    double expected_sigma_T = compute_expected_sigma_T(TARGET_TEMPERATURE_K, n_dof);

    printf("\n");
    printf("    === V.9: Temperature Fluctuation ===\n");
    printf("    Expected sigma_T (theory): %.2f K\n", expected_sigma_T);
    printf("    Acceptance: sigma_T < %.1f K\n", TEMPERATURE_FLUCTUATION_MAX_K);
    printf("\n");

    char xvg_path[512];
    snprintf(xvg_path, sizeof(xvg_path), "%s/nve_500ps/temperature.xvg", B5_DATA_DIR);

    FILE *f = fopen(xvg_path, "r");
    if (!f) {
        printf("    Temperature data not found: %s\n", xvg_path);
        printf("    SKIP: Temperature data required\n");
        return -1;
    }
    fclose(f);

    printf("    Temperature data found - analysis not yet implemented\n");
    printf("    SKIP: Full analysis pending implementation\n");
    return -1;
#endif
}

/* ===========================================================================
 * MAIN TEST RUNNER
 * ===========================================================================*/

static void print_usage(const char *prog)
{
    printf("Usage: %s [--unit | --integration | --all]\n", prog);
    printf("\n");
    printf("Options:\n");
    printf("  --unit         Run unit tests only (default)\n");
    printf("  --integration  Run integration tests only\n");
    printf("  --all          Run all tests\n");
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

    printf("=== US-043 Temperature Stability Tests (V.8-V.9) ===\n\n");

    int n_dof = compute_dof_periodic(B5_N_ATOMS_TOTAL);
    double expected_sigma_T = compute_expected_sigma_T(TARGET_TEMPERATURE_K, n_dof);

    printf("Acceptance Criteria (from T-US-043-V.8, V.9):\n");
    printf("  V.8: Temperature mean = %.0f +/- %.0f K\n",
           TARGET_TEMPERATURE_K, TEMPERATURE_MEAN_TOL_K);
    printf("  V.9: Temperature fluctuation < %.0f K\n", TEMPERATURE_FLUCTUATION_MAX_K);
    printf("\n");

    printf("B5 System Parameters:\n");
    printf("  Total atoms: %d\n", B5_N_ATOMS_TOTAL);
    printf("  Degrees of freedom: %d (periodic)\n", n_dof);
    printf("  Expected sigma_T: %.2f K (from T*sqrt(2/N_dof))\n", expected_sigma_T);
    printf("\n");

    if (run_unit) {
        printf("=== Unit Tests (Mathematical Identities) ===\n\n");

        printf("Mean Calculation:\n");
        RUN_TEST(test_temperature_mean_calculation_T_US_043_N_1);
        printf("\n");

        printf("Standard Deviation:\n");
        RUN_TEST(test_temperature_stddev_calculation_T_US_043_N_2);
        printf("\n");

        printf("Degrees of Freedom (T-US-043-V.9):\n");
        RUN_TEST(test_temperature_dof_periodic_T_US_043_V_9);
        printf("\n");

        printf("Fluctuation Formula (T-US-043-V.9):\n");
        RUN_TEST(test_temperature_fluctuation_formula_T_US_043_V_9);
        RUN_TEST(test_temperature_fluctuation_scaling_T_US_043_V_9);
        printf("\n");

        printf("Acceptance Bounds Validation:\n");
        RUN_TEST(test_temperature_bounds_sensible_T_US_043_V_8);
        printf("\n");
    }

    if (run_integration) {
        printf("=== Integration Tests (V.8-V.9 Trajectory Analysis) ===\n\n");

#ifndef GRODFTB_LONG_TESTS
        printf("NOTE: Integration tests require GRODFTB_LONG_TESTS=ON\n");
        printf("      Rebuild with: cmake -DGRODFTB_LONG_TESTS=ON ..\n\n");
#endif

        printf("V.8 - Temperature Mean:\n");
        RUN_TEST(test_temperature_mean_V_US_043_V_8);
        printf("\n");

        printf("V.9 - Temperature Fluctuation:\n");
        RUN_TEST(test_temperature_fluctuation_V_US_043_V_9);
        printf("\n");
    }

    printf("=== Summary ===\n");
    printf("Tests run:     %d\n", tests_run);
    printf("Tests passed:  %d\n", tests_passed);
    printf("Tests failed:  %d\n", tests_failed);
    printf("Tests skipped: %d\n", tests_skipped);

    if (tests_failed > 0) {
        printf("\nFAILED: %d test(s) failed\n", tests_failed);
        return 1;
    }

    return 0;
}
