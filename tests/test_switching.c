/*
 * SDD:specs.md:S8.1 -- Switching function mathematical identity tests
 * US-039: Switching Function for Cutoff Embedding
 *
 * ThDD:T-US-039-V.1 -- Switching function algebraic properties
 * ThDD:T-US-039-V.2 -- First derivative factorization
 * ThDD:T-US-039-V.3 -- Second derivative factorization
 *
 * These tests verify that the quintic switching function:
 *   S(u) = 1 - 10u^3 + 15u^4 - 6u^5
 * satisfies exact mathematical identities to machine precision.
 *
 * IMPORTANT: These are mathematical identity tests. The expected values
 * are exact algebraic results (S(0)=1, S(1)=0, S(0.5)=0.5, etc.) -- they
 * are NOT fabricated reference data. The tolerance of 10^-15 reflects
 * machine precision for double-precision floating-point arithmetic.
 *
 * Per CLAUDE.md Golden Rule: The expected values here are mathematical
 * facts, not empirical data requiring provenance from calculations.
 *
 * LGPL-3.0-or-later
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* ---------------------------------------------------------------------------
 * Switching function header (TO BE IMPLEMENTED)
 *
 * These functions do not exist yet -- this test is written FIRST per TDD.
 * The test will fail to link until the implementation is created.
 *
 * Expected API (SDD:specs.md:S8.1):
 *   double grodftb_switch_func(double u);
 *   double grodftb_switch_deriv(double u);
 *   double grodftb_switch_deriv2(double u);
 * --------------------------------------------------------------------------- */
#include "grodftb/switching.h"

/* ---------------------------------------------------------------------------
 * Test framework
 * --------------------------------------------------------------------------- */
static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

static void pass(const char *name)
{
    tests_run++;
    tests_passed++;
    printf("  PASS: %s\n", name);
}

static void fail(const char *name, const char *reason)
{
    tests_run++;
    tests_failed++;
    fprintf(stderr, "  FAIL: %s -- %s\n", name, reason);
}

/* ---------------------------------------------------------------------------
 * ThDD:T-US-039-V.1 -- Tolerance for mathematical identities
 *
 * For double precision (IEEE 754), machine epsilon is approximately 2.2e-16.
 * We use 10^-15 as the tolerance to allow for accumulated rounding error
 * in polynomial evaluation (5 terms, ~5 FLOPs).
 * --------------------------------------------------------------------------- */
#define IDENTITY_TOL 1e-15

/* ---------------------------------------------------------------------------
 * Test 1: test_switching_function_identities_T_US_039_V_1
 *
 * ThDD:T-US-039-V.1 -- Switching function endpoint and midpoint identities
 *
 * Mathematical facts (algebraic, not empirical):
 *   S(0) = 1      (full weight at inner boundary)
 *   S(1) = 0      (zero weight at outer boundary)
 *   S(0.5) = 0.5  (symmetry point)
 *
 * Derivation of S(0.5) = 0.5:
 *   S(0.5) = 1 - 10*(1/8) + 15*(1/16) - 6*(1/32)
 *          = 1 - 10/8 + 15/16 - 6/32
 *          = 32/32 - 40/32 + 30/32 - 6/32
 *          = 16/32 = 0.5
 * --------------------------------------------------------------------------- */
static void test_switching_function_identities_T_US_039_V_1(void)
{
    const char *test_name = "test_switching_function_identities_T_US_039_V_1";

    /* S(0) = 1 exactly */
    double s0 = grodftb_switch_func(0.0);
    double err_s0 = fabs(s0 - 1.0);
    if (err_s0 > IDENTITY_TOL) {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "S(0) = %.17e, expected 1.0, error = %.2e > %.2e",
                 s0, err_s0, IDENTITY_TOL);
        fail(test_name, msg);
        return;
    }

    /* S(1) = 0 exactly */
    double s1 = grodftb_switch_func(1.0);
    double err_s1 = fabs(s1 - 0.0);
    if (err_s1 > IDENTITY_TOL) {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "S(1) = %.17e, expected 0.0, error = %.2e > %.2e",
                 s1, err_s1, IDENTITY_TOL);
        fail(test_name, msg);
        return;
    }

    /* S(0.5) = 0.5 exactly (algebraic symmetry) */
    double s05 = grodftb_switch_func(0.5);
    double err_s05 = fabs(s05 - 0.5);
    if (err_s05 > IDENTITY_TOL) {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "S(0.5) = %.17e, expected 0.5, error = %.2e > %.2e",
                 s05, err_s05, IDENTITY_TOL);
        fail(test_name, msg);
        return;
    }

    pass(test_name);
}

/* ---------------------------------------------------------------------------
 * Test 2: test_switching_first_deriv_identities_T_US_039_V_1
 *
 * ThDD:T-US-039-V.1 -- First derivative boundary conditions
 *
 * Mathematical facts:
 *   S'(0) = 0  (continuous first derivative at inner boundary)
 *   S'(1) = 0  (continuous first derivative at outer boundary)
 *
 * This ensures force continuity when crossing the switching region boundary.
 * --------------------------------------------------------------------------- */
static void test_switching_first_deriv_identities_T_US_039_V_1(void)
{
    const char *test_name = "test_switching_first_deriv_identities_T_US_039_V_1";

    /* S'(0) = 0 exactly */
    double sp0 = grodftb_switch_deriv(0.0);
    double err_sp0 = fabs(sp0);
    if (err_sp0 > IDENTITY_TOL) {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "S'(0) = %.17e, expected 0.0, error = %.2e > %.2e",
                 sp0, err_sp0, IDENTITY_TOL);
        fail(test_name, msg);
        return;
    }

    /* S'(1) = 0 exactly */
    double sp1 = grodftb_switch_deriv(1.0);
    double err_sp1 = fabs(sp1);
    if (err_sp1 > IDENTITY_TOL) {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "S'(1) = %.17e, expected 0.0, error = %.2e > %.2e",
                 sp1, err_sp1, IDENTITY_TOL);
        fail(test_name, msg);
        return;
    }

    pass(test_name);
}

/* ---------------------------------------------------------------------------
 * Test 3: test_switching_second_deriv_identities_T_US_039_V_1
 *
 * ThDD:T-US-039-V.1 -- Second derivative boundary conditions
 *
 * Mathematical facts:
 *   S''(0) = 0  (continuous second derivative at inner boundary)
 *   S''(1) = 0  (continuous second derivative at outer boundary)
 *
 * This ensures smooth force variation (no kinks) at switching boundaries.
 * --------------------------------------------------------------------------- */
static void test_switching_second_deriv_identities_T_US_039_V_1(void)
{
    const char *test_name = "test_switching_second_deriv_identities_T_US_039_V_1";

    /* S''(0) = 0 exactly */
    double spp0 = grodftb_switch_deriv2(0.0);
    double err_spp0 = fabs(spp0);
    if (err_spp0 > IDENTITY_TOL) {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "S''(0) = %.17e, expected 0.0, error = %.2e > %.2e",
                 spp0, err_spp0, IDENTITY_TOL);
        fail(test_name, msg);
        return;
    }

    /* S''(1) = 0 exactly */
    double spp1 = grodftb_switch_deriv2(1.0);
    double err_spp1 = fabs(spp1);
    if (err_spp1 > IDENTITY_TOL) {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "S''(1) = %.17e, expected 0.0, error = %.2e > %.2e",
                 spp1, err_spp1, IDENTITY_TOL);
        fail(test_name, msg);
        return;
    }

    pass(test_name);
}

/* ---------------------------------------------------------------------------
 * Test 4: test_switching_derivative_factorization_T_US_039_V_2
 *
 * ThDD:T-US-039-V.2 -- First derivative factorization identity
 *
 * Mathematical identity:
 *   S'(u) = -30 u^2 (1-u)^2
 *
 * This factored form:
 * 1. Guarantees S'(0) = S'(1) = 0 algebraically
 * 2. Provides numerically stable computation
 * 3. Enables efficient force computation
 *
 * Derivation:
 *   S(u) = 1 - 10u^3 + 15u^4 - 6u^5
 *   S'(u) = -30u^2 + 60u^3 - 30u^4
 *        = -30u^2 (1 - 2u + u^2)
 *        = -30u^2 (1-u)^2
 *
 * Test: Compare grodftb_switch_deriv() against factored form at 100 points.
 * --------------------------------------------------------------------------- */
static void test_switching_derivative_factorization_T_US_039_V_2(void)
{
    const char *test_name = "test_switching_derivative_factorization_T_US_039_V_2";

    /* Test at 101 uniformly spaced points in [0, 1] */
    const int n_points = 101;
    const double factored_tol = 1e-14;  /* Slightly relaxed for accumulated error */

    for (int i = 0; i < n_points; i++) {
        double u = (double)i / (double)(n_points - 1);

        /* Factored form: S'(u) = -30 u^2 (1-u)^2 */
        double one_minus_u = 1.0 - u;
        double factored = -30.0 * u * u * one_minus_u * one_minus_u;

        /* Implementation value */
        double impl = grodftb_switch_deriv(u);

        /* Check agreement */
        double abs_err = fabs(impl - factored);
        double rel_err = (fabs(factored) > 1e-15)
                         ? abs_err / fabs(factored)
                         : abs_err;

        if (abs_err > factored_tol && rel_err > factored_tol) {
            char msg[256];
            snprintf(msg, sizeof(msg),
                     "At u=%.4f: impl=%.17e, factored=%.17e, "
                     "abs_err=%.2e, rel_err=%.2e",
                     u, impl, factored, abs_err, rel_err);
            fail(test_name, msg);
            return;
        }
    }

    pass(test_name);
}

/* ---------------------------------------------------------------------------
 * Test 5: test_switching_derivative_fd_T_US_039_V_2
 *
 * ThDD:T-US-039-V.2 -- Validate derivative against finite differences
 *
 * Central difference: S'(u) approx [S(u+d) - S(u-d)] / (2d)
 * with d = 10^-6.
 *
 * This verifies the implementation without relying on the factored form.
 * --------------------------------------------------------------------------- */
static void test_switching_derivative_fd_T_US_039_V_2(void)
{
    const char *test_name = "test_switching_derivative_fd_T_US_039_V_2";

    const double delta = 1e-6;
    const double fd_tol = 1e-8;  /* FD truncation error O(delta^2) ~ 10^-12 */
    const int n_points = 21;     /* Test at 21 interior points */

    for (int i = 1; i < n_points - 1; i++) {
        double u = (double)i / (double)(n_points - 1);

        /* Skip boundary points where derivative is zero */
        if (u < delta || u > 1.0 - delta) continue;

        /* Finite difference approximation */
        double s_plus = grodftb_switch_func(u + delta);
        double s_minus = grodftb_switch_func(u - delta);
        double fd_deriv = (s_plus - s_minus) / (2.0 * delta);

        /* Analytic derivative */
        double analytic = grodftb_switch_deriv(u);

        /* Check agreement */
        double abs_err = fabs(analytic - fd_deriv);
        double rel_err = (fabs(analytic) > 1e-15)
                         ? abs_err / fabs(analytic)
                         : abs_err;

        if (rel_err > fd_tol && abs_err > fd_tol) {
            char msg[256];
            snprintf(msg, sizeof(msg),
                     "At u=%.4f: analytic=%.10e, FD=%.10e, rel_err=%.2e",
                     u, analytic, fd_deriv, rel_err);
            fail(test_name, msg);
            return;
        }
    }

    pass(test_name);
}

/* ---------------------------------------------------------------------------
 * Test 6: test_switching_second_deriv_factorization_T_US_039_V_3
 *
 * ThDD:T-US-039-V.3 -- Second derivative factorization identity
 *
 * Mathematical identity:
 *   S''(u) = -60 u (1-u) (1-2u)
 *
 * Properties:
 *   S''(0) = 0 (factor u = 0)
 *   S''(1) = 0 (factor (1-u) = 0)
 *   S''(0.5) = 0 (factor (1-2u) = 0, inflection point)
 *
 * Sign analysis for S''(u) = -60 u (1-u) (1-2u):
 *   - The factor "-60" is negative
 *   - u > 0 for u in (0, 1)
 *   - (1-u) > 0 for u in (0, 1)
 *   - (1-2u) > 0 for u in (0, 0.5), = 0 at u=0.5, < 0 for u in (0.5, 1)
 *
 * Therefore:
 *   S''(u) < 0 for u in (0, 0.5): slope S' is decreasing (becoming more negative)
 *   S''(u) > 0 for u in (0.5, 1): slope S' is increasing (becoming less negative)
 *
 * This matches the shape of S(u): the switching function has its steepest
 * downward slope at u = 0.5, and the slope is zero at both boundaries.
 * --------------------------------------------------------------------------- */
static void test_switching_second_deriv_factorization_T_US_039_V_3(void)
{
    const char *test_name = "test_switching_second_deriv_factorization_T_US_039_V_3";

    const double factored_tol = 1e-14;
    const int n_points = 101;

    for (int i = 0; i < n_points; i++) {
        double u = (double)i / (double)(n_points - 1);

        /* Factored form: S''(u) = -60 u (1-u) (1-2u) */
        double factored = -60.0 * u * (1.0 - u) * (1.0 - 2.0 * u);

        /* Implementation value */
        double impl = grodftb_switch_deriv2(u);

        /* Check agreement */
        double abs_err = fabs(impl - factored);

        if (abs_err > factored_tol) {
            char msg[256];
            snprintf(msg, sizeof(msg),
                     "At u=%.4f: impl=%.17e, factored=%.17e, abs_err=%.2e",
                     u, impl, factored, abs_err);
            fail(test_name, msg);
            return;
        }
    }

    /* Verify sign pattern:
     * S''(0.25) < 0: slope S' is still decreasing (becoming more negative)
     * S''(0.75) > 0: slope S' is increasing (becoming less negative)
     */
    double spp_025 = grodftb_switch_deriv2(0.25);
    double spp_075 = grodftb_switch_deriv2(0.75);

    if (spp_025 >= 0.0) {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "S''(0.25) = %.10e, expected < 0 (slope still decreasing)", spp_025);
        fail(test_name, msg);
        return;
    }
    if (spp_075 <= 0.0) {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "S''(0.75) = %.10e, expected > 0 (slope increasing)", spp_075);
        fail(test_name, msg);
        return;
    }

    /* Verify inflection at u = 0.5 */
    double spp_05 = grodftb_switch_deriv2(0.5);
    if (fabs(spp_05) > IDENTITY_TOL) {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "S''(0.5) = %.17e, expected 0.0 (inflection point)",
                 spp_05);
        fail(test_name, msg);
        return;
    }

    pass(test_name);
}

/* ---------------------------------------------------------------------------
 * Test 7: test_switching_monotonicity_T_US_039_V_1
 *
 * ThDD:T-US-039-V.1 -- Switching function monotonicity
 *
 * Mathematical property: S(u) is monotonically decreasing on [0, 1].
 * Since S'(u) = -30 u^2 (1-u)^2 <= 0 for all u in [0, 1], with equality
 * only at the endpoints.
 * --------------------------------------------------------------------------- */
static void test_switching_monotonicity_T_US_039_V_1(void)
{
    const char *test_name = "test_switching_monotonicity_T_US_039_V_1";

    const int n_points = 1001;
    double prev_s = grodftb_switch_func(0.0);

    for (int i = 1; i < n_points; i++) {
        double u = (double)i / (double)(n_points - 1);
        double s = grodftb_switch_func(u);

        if (s > prev_s + IDENTITY_TOL) {
            char msg[256];
            snprintf(msg, sizeof(msg),
                     "Monotonicity violated: S(%.6f)=%.10e > S(%.6f)=%.10e",
                     u, s, (double)(i-1)/(double)(n_points-1), prev_s);
            fail(test_name, msg);
            return;
        }

        prev_s = s;
    }

    pass(test_name);
}

/* ---------------------------------------------------------------------------
 * Test 8: test_switching_range_T_US_039_V_1
 *
 * ThDD:T-US-039-V.1 -- Switching function stays in [0, 1]
 *
 * Mathematical property: For u in [0, 1], we have S(u) in [0, 1].
 * --------------------------------------------------------------------------- */
static void test_switching_range_T_US_039_V_1(void)
{
    const char *test_name = "test_switching_range_T_US_039_V_1";

    const int n_points = 1001;

    for (int i = 0; i < n_points; i++) {
        double u = (double)i / (double)(n_points - 1);
        double s = grodftb_switch_func(u);

        if (s < -IDENTITY_TOL || s > 1.0 + IDENTITY_TOL) {
            char msg[256];
            snprintf(msg, sizeof(msg),
                     "Range violation: S(%.6f) = %.17e not in [0, 1]",
                     u, s);
            fail(test_name, msg);
            return;
        }
    }

    pass(test_name);
}

/* ---------------------------------------------------------------------------
 * main
 * --------------------------------------------------------------------------- */
int main(void)
{
    printf("=== US-039 Switching Function Identity Tests ===\n\n");
    printf("ThDD:T-US-039-V.1 -- Algebraic properties of quintic switch\n");
    printf("ThDD:T-US-039-V.2 -- First derivative factorization\n");
    printf("ThDD:T-US-039-V.3 -- Second derivative factorization\n\n");

    printf("Note: These test mathematical identities, not empirical data.\n");
    printf("Expected values (S(0)=1, S(1)=0, S(0.5)=0.5, etc.) are exact.\n\n");

    /* Run tests */
    test_switching_function_identities_T_US_039_V_1();
    test_switching_first_deriv_identities_T_US_039_V_1();
    test_switching_second_deriv_identities_T_US_039_V_1();
    test_switching_derivative_factorization_T_US_039_V_2();
    test_switching_derivative_fd_T_US_039_V_2();
    test_switching_second_deriv_factorization_T_US_039_V_3();
    test_switching_monotonicity_T_US_039_V_1();
    test_switching_range_T_US_039_V_1();

    /* Summary */
    printf("\n=== Summary ===\n");
    printf("Tests run: %d, passed: %d, failed: %d\n",
           tests_run, tests_passed, tests_failed);

    return (tests_failed > 0) ? 1 : 0;
}
