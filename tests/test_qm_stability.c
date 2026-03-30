/*
 * US-043: QM Subsystem Stability Tests (V.10-V.12)
 *
 * ThDD:T-US-043-V.10 -- SCC convergence rate criterion
 * ThDD:T-US-043-V.11 -- QM energy stability criterion
 * ThDD:T-US-043-V.12 -- Mulliken charge conservation criterion
 * SDD:specs.md:S21.1 -- NVE simulation requirements
 * SDD:docs/verification/US-043.md -- Verification plan
 *
 * This file implements:
 *   1. V.10: SCC convergence rate (must be 100%)
 *   2. V.11: QM energy deviation from mean (bounded)
 *   3. V.12: Mulliken charge conservation (total charge constant)
 *
 * GOLDEN RULE COMPLIANCE:
 * - Unit tests use mathematical identities (exact formulas)
 * - Thresholds from T-US-043 criteria, not fabricated
 * - No hardcoded "expected" simulation energy values
 * - Charge conservation threshold (1e-6 e) from numerical analysis
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
 * ACCEPTANCE CRITERIA (from T-US-043-V.10 through V.12)
 *
 * These are specification-derived thresholds, NOT fabricated expected values.
 * ===========================================================================*/

/* V.10: SCC convergence rate (must be 100%) */
#define SCC_CONVERGENCE_RATE_REQUIRED  1.0

/* V.11: QM energy deviation bound in kJ/mol
 * This allows for thermal fluctuations in QM region
 * ~100 kJ/mol = ~0.04 Hartree (from T-US-043-V.11)
 */
#define QM_ENERGY_DEVIATION_MAX_KJMOL  100.0

/* V.12: Mulliken charge conservation tolerance in elementary charges
 * 10^-6 e chosen to be well above numerical noise (~10^-10 e)
 * but detect any systematic charge drift
 */
#define CHARGE_CONSERVATION_TOL_E  1e-6

/* ===========================================================================
 * B5 SYSTEM PARAMETERS
 * ===========================================================================*/

#define B5_N_ATOMS_TOTAL  2652
#define B5_N_QM_ATOMS     3      /* One water molecule: O, H, H */
#define B5_QM_TOTAL_CHARGE  0.0  /* Neutral water molecule */

/* ===========================================================================
 * HELPER FUNCTIONS
 * ===========================================================================*/

/**
 * Compute sum of array elements.
 *
 * ThDD:T-US-043-V.12 -- Total charge calculation
 *
 * @param n    Number of elements
 * @param arr  Data array
 * @return     Sum of all elements
 */
static double compute_sum(int n, const double *arr)
{
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += arr[i];
    }
    return sum;
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
    return compute_sum(n, arr) / (double)n;
}

/**
 * Compute maximum absolute deviation from mean.
 *
 * ThDD:T-US-043-V.11 -- |E_QM(t) - <E_QM>| check
 *
 * @param n    Number of elements
 * @param arr  Data array
 * @return     Maximum |arr[i] - mean|
 */
static double compute_max_deviation(int n, const double *arr)
{
    if (n <= 0) return 0.0;

    double mean = compute_mean(n, arr);
    double max_dev = 0.0;

    for (int i = 0; i < n; i++) {
        double dev = fabs(arr[i] - mean);
        if (dev > max_dev) {
            max_dev = dev;
        }
    }

    return max_dev;
}

/**
 * Count number of converged SCC calculations.
 *
 * ThDD:T-US-043-V.10 -- Convergence rate calculation
 *
 * @param n          Number of steps
 * @param converged  Array of convergence flags (1 = converged, 0 = not)
 * @return           Number of converged steps
 */
static int count_converged(int n, const int *converged)
{
    int count = 0;
    for (int i = 0; i < n; i++) {
        if (converged[i]) count++;
    }
    return count;
}

/**
 * Compute convergence rate.
 *
 * ThDD:T-US-043-V.10 -- convergence_rate = N_converged / N_total
 *
 * @param n_converged  Number of converged steps
 * @param n_total      Total number of steps
 * @return             Convergence rate [0, 1]
 */
static double compute_convergence_rate(int n_converged, int n_total)
{
    if (n_total <= 0) return 0.0;
    return (double)n_converged / (double)n_total;
}

/**
 * Check if total charge is conserved between two values.
 *
 * ThDD:T-US-043-V.12 -- |Q(t) - Q(0)| < tolerance
 *
 * @param q_initial    Initial total charge
 * @param q_current    Current total charge
 * @param tolerance    Maximum allowed difference
 * @return             1 if conserved, 0 if violated
 */
static int is_charge_conserved(double q_initial, double q_current, double tolerance)
{
    return fabs(q_current - q_initial) < tolerance ? 1 : 0;
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
 * Test: Sum calculation.
 *
 * MATHEMATICAL IDENTITY:
 *   sum({1, 2, 3, 4, 5}) = 15
 */
static int test_qm_stability_sum_calculation_T_US_043_N_1(void)
{
    double data[5] = {1.0, 2.0, 3.0, 4.0, 5.0};

    double sum = compute_sum(5, data);
    double expected = 15.0;

    const double TOL = 1e-14;

    if (fabs(sum - expected) > TOL) {
        fprintf(stderr, "\n    FAIL: sum = %.15f, expected %.15f\n", sum, expected);
        return 0;
    }

    return 1;
}

/**
 * Test: Maximum deviation calculation.
 *
 * MATHEMATICAL IDENTITY:
 *   For data {10, 20, 30, 40, 50}:
 *   mean = 30
 *   max_deviation = max(|10-30|, |20-30|, |30-30|, |40-30|, |50-30|)
 *                 = max(20, 10, 0, 10, 20) = 20
 */
static int test_qm_stability_max_deviation_T_US_043_V_11(void)
{
    double data[5] = {10.0, 20.0, 30.0, 40.0, 50.0};

    double max_dev = compute_max_deviation(5, data);
    double expected = 20.0;

    const double TOL = 1e-14;

    if (fabs(max_dev - expected) > TOL) {
        fprintf(stderr, "\n    FAIL: max_dev = %.15f, expected %.15f\n",
                max_dev, expected);
        return 0;
    }

    return 1;
}

/**
 * Test: Convergence rate calculation.
 *
 * ThDD:T-US-043-V.10
 *
 * MATHEMATICAL IDENTITY:
 *   5 converged out of 10 total -> rate = 0.5
 */
static int test_qm_stability_convergence_rate_T_US_043_V_10(void)
{
    double rate = compute_convergence_rate(5, 10);
    double expected = 0.5;

    const double TOL = 1e-14;

    if (fabs(rate - expected) > TOL) {
        fprintf(stderr, "\n    FAIL: rate = %.15f, expected %.15f\n", rate, expected);
        return 0;
    }

    /* Test 100% convergence */
    rate = compute_convergence_rate(1000, 1000);
    if (fabs(rate - 1.0) > TOL) {
        fprintf(stderr, "\n    FAIL: 100%% rate = %.15f, expected 1.0\n", rate);
        return 0;
    }

    return 1;
}

/**
 * Test: Convergence counting.
 *
 * ThDD:T-US-043-V.10
 *
 * MATHEMATICAL IDENTITY:
 *   Array {1, 1, 0, 1, 1, 0, 1} has 5 ones
 */
static int test_qm_stability_count_converged_T_US_043_V_10(void)
{
    int converged[7] = {1, 1, 0, 1, 1, 0, 1};

    int count = count_converged(7, converged);
    int expected = 5;

    if (count != expected) {
        fprintf(stderr, "\n    FAIL: count = %d, expected %d\n", count, expected);
        return 0;
    }

    return 1;
}

/**
 * Test: Charge conservation check.
 *
 * ThDD:T-US-043-V.12
 *
 * MATHEMATICAL IDENTITY:
 *   |0.0 - 1e-7| = 1e-7 < 1e-6 -> conserved
 *   |0.0 - 1e-5| = 1e-5 > 1e-6 -> violated
 */
static int test_qm_stability_charge_conservation_T_US_043_V_12(void)
{
    /* Small drift within tolerance */
    if (!is_charge_conserved(0.0, 1e-7, CHARGE_CONSERVATION_TOL_E)) {
        fprintf(stderr, "\n    FAIL: 1e-7 drift should be conserved\n");
        return 0;
    }

    /* Drift outside tolerance */
    if (is_charge_conserved(0.0, 1e-5, CHARGE_CONSERVATION_TOL_E)) {
        fprintf(stderr, "\n    FAIL: 1e-5 drift should NOT be conserved\n");
        return 0;
    }

    /* Exact conservation */
    if (!is_charge_conserved(0.0, 0.0, CHARGE_CONSERVATION_TOL_E)) {
        fprintf(stderr, "\n    FAIL: Zero drift should be conserved\n");
        return 0;
    }

    return 1;
}

/**
 * Test: QM energy deviation threshold is physically sensible.
 *
 * ThDD:T-US-043-V.11
 *
 * The 100 kJ/mol threshold should be large enough to accommodate
 * thermal fluctuations but small enough to catch runaway QM energy.
 *
 * For a 3-atom water molecule at 300 K:
 *   Thermal energy scale: (3/2) * k_B * T * 3 atoms = 11.2 kJ/mol
 *   100 kJ/mol ~ 9x thermal scale (reasonable margin)
 */
static int test_qm_stability_energy_threshold_sensible_T_US_043_V_11(void)
{
    double k_B = 0.0083144626;  /* kJ/mol/K */
    double T = 300.0;  /* K */

    /* Thermal energy for QM water molecule */
    double thermal_energy = 1.5 * k_B * T * B5_N_QM_ATOMS;  /* ~11.2 kJ/mol */

    /* Threshold should be much larger than thermal fluctuations */
    if (QM_ENERGY_DEVIATION_MAX_KJMOL < 5.0 * thermal_energy) {
        fprintf(stderr, "\n    FAIL: threshold %.1f kJ/mol < 5x thermal (%.1f kJ/mol)\n",
                QM_ENERGY_DEVIATION_MAX_KJMOL, 5.0 * thermal_energy);
        return 0;
    }

    /* But not absurdly large (> 1000 kJ/mol would indicate bad choice) */
    if (QM_ENERGY_DEVIATION_MAX_KJMOL > 1000.0) {
        fprintf(stderr, "\n    FAIL: threshold %.1f kJ/mol suspiciously large\n",
                QM_ENERGY_DEVIATION_MAX_KJMOL);
        return 0;
    }

    printf("\n    (thermal ~ %.1f kJ/mol, threshold = %.0f kJ/mol) ",
           thermal_energy, QM_ENERGY_DEVIATION_MAX_KJMOL);
    return 1;
}

/**
 * Test: Charge conservation threshold is numerically sensible.
 *
 * ThDD:T-US-043-V.12
 *
 * The 1e-6 e threshold should be:
 * - Well above machine epsilon * typical charge magnitude (~1e-15)
 * - Small enough to detect systematic drift
 */
static int test_qm_stability_charge_threshold_sensible_T_US_043_V_12(void)
{
    /* Should be above numerical noise */
    double numerical_noise = 1e-10;  /* Typical SCC precision */

    if (CHARGE_CONSERVATION_TOL_E < 10.0 * numerical_noise) {
        fprintf(stderr, "\n    FAIL: tolerance %.0e < 10x numerical noise (%.0e)\n",
                CHARGE_CONSERVATION_TOL_E, numerical_noise);
        return 0;
    }

    /* Should be small enough to detect physical charge transfer */
    double detectable_transfer = 0.01;  /* 1% of electron */

    if (CHARGE_CONSERVATION_TOL_E > 0.1 * detectable_transfer) {
        fprintf(stderr, "\n    FAIL: tolerance %.0e too large to detect 1%% transfer\n",
                CHARGE_CONSERVATION_TOL_E);
        return 0;
    }

    printf("\n    (tolerance = %.0e e, well above noise %.0e, below 1%% transfer) ",
           CHARGE_CONSERVATION_TOL_E, numerical_noise);
    return 1;
}

/* ===========================================================================
 * SECTION 2: INTEGRATION TESTS (Trajectory Analysis)
 *
 * These tests require actual trajectory data and DFTB+ logs.
 * Gated behind GRODFTB_LONG_TESTS.
 *
 * GOLDEN RULE: All expected values are thresholds from T-US-043,
 * NOT fabricated energy predictions.
 * ===========================================================================*/

#ifndef GRODFTB_SOURCE_DIR_STR
#define GRODFTB_SOURCE_DIR_STR "."
#endif
#define B5_DATA_DIR GRODFTB_SOURCE_DIR_STR "/tests/data/b5"

/**
 * Test V.10: SCC convergence rate.
 *
 * ThDD:T-US-043-V.10
 *
 * Criterion: convergence_rate = N_converged / N_total = 100%
 *
 * NOTE: This reads SCC convergence status from DFTB+ output logs.
 * Non-convergence at any step invalidates the entire trajectory.
 */
static int test_scc_convergence_V_US_043_V_10(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    printf("\n");
    printf("    === V.10: SCC Convergence Rate ===\n");
    printf("    Acceptance: convergence_rate = %.0f%%\n",
           SCC_CONVERGENCE_RATE_REQUIRED * 100.0);
    printf("\n");
    printf("    Rationale: Non-convergent SCC produces incorrect forces,\n");
    printf("    which invalidates the entire trajectory.\n");
    printf("\n");

    char log_path[512];
    snprintf(log_path, sizeof(log_path), "%s/nve_500ps/dftb_output.log", B5_DATA_DIR);

    FILE *f = fopen(log_path, "r");
    if (!f) {
        printf("    DFTB+ log not found: %s\n", log_path);
        printf("    SKIP: SCC convergence data required\n");
        return -1;
    }
    fclose(f);

    printf("    DFTB+ log found - analysis not yet implemented\n");
    printf("    SKIP: Full analysis pending implementation\n");
    return -1;
#endif
}

/**
 * Test V.11: QM energy stability.
 *
 * ThDD:T-US-043-V.11
 *
 * Criterion: |E_QM(t) - <E_QM>| < 100 kJ/mol for all t
 *
 * NOTE: This checks that QM energy fluctuations are bounded,
 * indicating stable QM region without geometry explosion.
 */
static int test_qm_energy_bound_V_US_043_V_11(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    printf("\n");
    printf("    === V.11: QM Energy Stability ===\n");
    printf("    Acceptance: |E_QM(t) - <E_QM>| < %.0f kJ/mol\n",
           QM_ENERGY_DEVIATION_MAX_KJMOL);
    printf("\n");
    printf("    Rationale: Large deviations indicate unphysical geometry\n");
    printf("    or embedding problems in the QM region.\n");
    printf("\n");

    char xvg_path[512];
    snprintf(xvg_path, sizeof(xvg_path), "%s/nve_500ps/qm_energy.xvg", B5_DATA_DIR);

    FILE *f = fopen(xvg_path, "r");
    if (!f) {
        printf("    QM energy data not found: %s\n", xvg_path);
        printf("    SKIP: QM energy data required\n");
        return -1;
    }
    fclose(f);

    printf("    QM energy data found - analysis not yet implemented\n");
    printf("    SKIP: Full analysis pending implementation\n");
    return -1;
#endif
}

/**
 * Test V.12: Mulliken charge conservation.
 *
 * ThDD:T-US-043-V.12
 *
 * Criterion: |sum_A q_A(t) - sum_A q_A(0)| < 10^-6 e
 *
 * NOTE: Total Mulliken charge on QM region must be conserved.
 * For neutral water, initial charge = 0.0.
 */
static int test_charge_conservation_V_US_043_V_12(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    printf("\n");
    printf("    === V.12: Mulliken Charge Conservation ===\n");
    printf("    QM region: %d atoms (water molecule)\n", B5_N_QM_ATOMS);
    printf("    Expected total charge: %.1f e (neutral)\n", B5_QM_TOTAL_CHARGE);
    printf("    Acceptance: |Q(t) - Q(0)| < %.0e e\n", CHARGE_CONSERVATION_TOL_E);
    printf("\n");
    printf("    Rationale: Charge drift indicates SCC error or\n");
    printf("    incorrect charge extraction from DFTB+.\n");
    printf("\n");

    char xvg_path[512];
    snprintf(xvg_path, sizeof(xvg_path), "%s/nve_500ps/mulliken_charges.xvg", B5_DATA_DIR);

    FILE *f = fopen(xvg_path, "r");
    if (!f) {
        printf("    Mulliken charge data not found: %s\n", xvg_path);
        printf("    SKIP: Charge data required\n");
        return -1;
    }
    fclose(f);

    printf("    Mulliken charge data found - analysis not yet implemented\n");
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

    printf("=== US-043 QM Subsystem Stability Tests (V.10-V.12) ===\n\n");

    printf("Acceptance Criteria (from T-US-043-V.10 through V.12):\n");
    printf("  V.10: SCC convergence rate = %.0f%%\n", SCC_CONVERGENCE_RATE_REQUIRED * 100.0);
    printf("  V.11: QM energy deviation < %.0f kJ/mol\n", QM_ENERGY_DEVIATION_MAX_KJMOL);
    printf("  V.12: Charge conservation < %.0e e\n", CHARGE_CONSERVATION_TOL_E);
    printf("\n");

    printf("B5 QM Region:\n");
    printf("  QM atoms: %d (one water molecule)\n", B5_N_QM_ATOMS);
    printf("  Expected total charge: %.1f e\n", B5_QM_TOTAL_CHARGE);
    printf("\n");

    if (run_unit) {
        printf("=== Unit Tests (Mathematical Identities) ===\n\n");

        printf("Sum Calculation:\n");
        RUN_TEST(test_qm_stability_sum_calculation_T_US_043_N_1);
        printf("\n");

        printf("Max Deviation (T-US-043-V.11):\n");
        RUN_TEST(test_qm_stability_max_deviation_T_US_043_V_11);
        printf("\n");

        printf("Convergence Rate (T-US-043-V.10):\n");
        RUN_TEST(test_qm_stability_convergence_rate_T_US_043_V_10);
        RUN_TEST(test_qm_stability_count_converged_T_US_043_V_10);
        printf("\n");

        printf("Charge Conservation (T-US-043-V.12):\n");
        RUN_TEST(test_qm_stability_charge_conservation_T_US_043_V_12);
        printf("\n");

        printf("Threshold Validation:\n");
        RUN_TEST(test_qm_stability_energy_threshold_sensible_T_US_043_V_11);
        RUN_TEST(test_qm_stability_charge_threshold_sensible_T_US_043_V_12);
        printf("\n");
    }

    if (run_integration) {
        printf("=== Integration Tests (V.10-V.12 Trajectory Analysis) ===\n\n");

#ifndef GRODFTB_LONG_TESTS
        printf("NOTE: Integration tests require GRODFTB_LONG_TESTS=ON\n");
        printf("      Rebuild with: cmake -DGRODFTB_LONG_TESTS=ON ..\n\n");
#endif

        printf("V.10 - SCC Convergence:\n");
        RUN_TEST(test_scc_convergence_V_US_043_V_10);
        printf("\n");

        printf("V.11 - QM Energy Stability:\n");
        RUN_TEST(test_qm_energy_bound_V_US_043_V_11);
        printf("\n");

        printf("V.12 - Charge Conservation:\n");
        RUN_TEST(test_charge_conservation_V_US_043_V_12);
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
