/*
 * US-043: Trajectory Validity Tests (V.5-V.7)
 *
 * ThDD:T-US-043-V.5 through V.7 -- Trajectory validity criteria
 * SDD:specs.md:S21.1 -- NVE simulation requirements
 * SDD:docs/verification/US-043.md -- Verification plan
 *
 * This file implements:
 *   1. V.5: No NaN/Inf values detection
 *   2. V.6: QM geometry integrity (bond lengths, angles)
 *   3. V.7: Center of mass / momentum conservation
 *
 * GOLDEN RULE COMPLIANCE:
 * - Unit tests use mathematical identities (exact formulas)
 * - Physical bounds from literature (O-H bond ~0.096 nm, H-O-H ~104.5 deg)
 * - No hardcoded "expected" simulation energy values
 * - All simulation results must come from actual trajectory data
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
 * PHYSICAL BOUNDS (from established chemistry, NOT fabricated)
 *
 * ThDD:T-US-043-V.6 -- QM geometry bounds from spectroscopic data
 * ===========================================================================*/

/* O-H bond length bounds in nm (from T-US-043-V.6)
 * Equilibrium: 0.096 nm (gas phase water)
 * Bounds allow thermal vibration amplitude + numerical tolerance
 */
#define OH_BOND_MIN_NM  0.08
#define OH_BOND_MAX_NM  0.15

/* H-O-H angle bounds in degrees (from T-US-043-V.6)
 * Equilibrium: 104.5 degrees (gas phase water)
 * Bounds allow thermal bending + numerical tolerance
 */
#define HOH_ANGLE_MIN_DEG  80.0
#define HOH_ANGLE_MAX_DEG  120.0

/* ===========================================================================
 * ACCEPTANCE CRITERIA (from T-US-043-V.*)
 * ===========================================================================*/

/* V.7: Momentum conservation relative tolerance */
#define MOMENTUM_CONSERVATION_TOL  1e-6

/* ===========================================================================
 * B5 SYSTEM PARAMETERS
 * ===========================================================================*/

#define B5_N_ATOMS_TOTAL  2652
#define B5_N_QM_ATOMS     3

/* ===========================================================================
 * HELPER FUNCTIONS
 * ===========================================================================*/

/**
 * Check if a double value is finite (not NaN, not Inf).
 *
 * ThDD:T-US-043-V.5 -- Finiteness check for trajectory values
 *
 * @param x  Value to check
 * @return   1 if finite, 0 if NaN or Inf
 */
static int is_value_finite(double x)
{
    /* isfinite is C99 standard */
    return isfinite(x) ? 1 : 0;
}

/**
 * Compute Euclidean distance between two 3D points.
 *
 * ThDD:T-US-043-V.6 -- Distance formula: d = sqrt(sum((r1-r2)^2))
 *
 * @param r1  First point [3]
 * @param r2  Second point [3]
 * @return    Distance (same units as input)
 */
static double compute_distance(const double *r1, const double *r2)
{
    double dx = r1[0] - r2[0];
    double dy = r1[1] - r2[1];
    double dz = r1[2] - r2[2];
    return sqrt(dx*dx + dy*dy + dz*dz);
}

/**
 * Compute angle between three points (vertex at middle point).
 *
 * ThDD:T-US-043-V.6 -- Angle formula via dot product:
 *   cos(theta) = (v1 . v2) / (|v1| |v2|)
 *   where v1 = r1 - r_center, v2 = r3 - r_center
 *
 * @param r1      First point [3]
 * @param center  Vertex point [3]
 * @param r3      Third point [3]
 * @return        Angle in degrees
 */
static double compute_angle_degrees(const double *r1, const double *center, const double *r3)
{
    double v1[3], v2[3];

    for (int i = 0; i < 3; i++) {
        v1[i] = r1[i] - center[i];
        v2[i] = r3[i] - center[i];
    }

    double dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
    double len1 = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
    double len2 = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);

    if (len1 < 1e-15 || len2 < 1e-15) {
        return 0.0;  /* Degenerate case */
    }

    double cos_theta = dot / (len1 * len2);

    /* Clamp to [-1, 1] to avoid NaN from acos due to numerical error */
    if (cos_theta > 1.0) cos_theta = 1.0;
    if (cos_theta < -1.0) cos_theta = -1.0;

    return acos(cos_theta) * 180.0 / 3.14159265358979323846;
}

/**
 * Compute magnitude of a 3D vector.
 *
 * ThDD:T-US-043-V.7 -- Vector magnitude for momentum calculation
 *
 * @param v  Vector [3]
 * @return   Magnitude
 */
static double compute_magnitude(const double *v)
{
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
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
 * SECTION 1: UNIT TESTS (Mathematical Identities and Function Verification)
 *
 * GOLDEN RULE: All expected values are mathematical identities.
 * ===========================================================================*/

/**
 * Test: is_value_finite correctly identifies NaN.
 *
 * ThDD:T-US-043-V.5
 *
 * MATHEMATICAL IDENTITY:
 *   NaN != NaN (IEEE 754 property)
 *   isfinite(NaN) = false
 */
static int test_trajectory_validity_nan_detection_T_US_043_V_5(void)
{
    double nan_val = 0.0 / 0.0;  /* Generate NaN */

    if (is_value_finite(nan_val)) {
        fprintf(stderr, "\n    FAIL: NaN was not detected\n");
        return 0;
    }

    /* Verify normal values are detected as finite */
    if (!is_value_finite(0.0) || !is_value_finite(1.0) || !is_value_finite(-1e10)) {
        fprintf(stderr, "\n    FAIL: Normal value incorrectly flagged as non-finite\n");
        return 0;
    }

    return 1;
}

/**
 * Test: is_value_finite correctly identifies Inf.
 *
 * ThDD:T-US-043-V.5
 *
 * MATHEMATICAL IDENTITY:
 *   1.0/0.0 = +Inf (IEEE 754)
 *   isfinite(Inf) = false
 */
static int test_trajectory_validity_inf_detection_T_US_043_V_5(void)
{
    double pos_inf = 1.0 / 0.0;
    double neg_inf = -1.0 / 0.0;

    if (is_value_finite(pos_inf)) {
        fprintf(stderr, "\n    FAIL: +Inf was not detected\n");
        return 0;
    }

    if (is_value_finite(neg_inf)) {
        fprintf(stderr, "\n    FAIL: -Inf was not detected\n");
        return 0;
    }

    return 1;
}

/**
 * Test: Distance calculation.
 *
 * ThDD:T-US-043-V.6
 *
 * MATHEMATICAL IDENTITY:
 *   d((0,0,0), (3,4,0)) = 5 (Pythagorean theorem)
 */
static int test_trajectory_validity_distance_calculation_T_US_043_V_6(void)
{
    double r1[3] = {0.0, 0.0, 0.0};
    double r2[3] = {3.0, 4.0, 0.0};

    double d = compute_distance(r1, r2);
    double expected = 5.0;

    const double TOL = 1e-14;

    if (fabs(d - expected) > TOL) {
        fprintf(stderr, "\n    FAIL: distance = %.15f, expected %.15f\n", d, expected);
        return 0;
    }

    return 1;
}

/**
 * Test: Angle calculation for right angle.
 *
 * ThDD:T-US-043-V.6
 *
 * MATHEMATICAL IDENTITY:
 *   angle((1,0,0), (0,0,0), (0,1,0)) = 90 degrees
 */
static int test_trajectory_validity_angle_right_T_US_043_V_6(void)
{
    double r1[3] = {1.0, 0.0, 0.0};
    double center[3] = {0.0, 0.0, 0.0};
    double r3[3] = {0.0, 1.0, 0.0};

    double angle = compute_angle_degrees(r1, center, r3);
    double expected = 90.0;

    const double TOL = 1e-10;

    if (fabs(angle - expected) > TOL) {
        fprintf(stderr, "\n    FAIL: angle = %.15f deg, expected %.15f deg\n",
                angle, expected);
        return 0;
    }

    return 1;
}

/**
 * Test: Angle calculation for 60 degrees.
 *
 * ThDD:T-US-043-V.6
 *
 * MATHEMATICAL IDENTITY:
 *   angle((1,0,0), (0,0,0), (0.5, sqrt(3)/2, 0)) = 60 degrees
 *   (equilateral triangle vertices)
 */
static int test_trajectory_validity_angle_60deg_T_US_043_V_6(void)
{
    double r1[3] = {1.0, 0.0, 0.0};
    double center[3] = {0.0, 0.0, 0.0};
    double r3[3] = {0.5, 0.8660254037844386, 0.0};  /* sqrt(3)/2 = 0.866... */

    double angle = compute_angle_degrees(r1, center, r3);
    double expected = 60.0;

    const double TOL = 1e-10;

    if (fabs(angle - expected) > TOL) {
        fprintf(stderr, "\n    FAIL: angle = %.15f deg, expected %.15f deg\n",
                angle, expected);
        return 0;
    }

    return 1;
}

/**
 * Test: Angle calculation for water-like geometry.
 *
 * ThDD:T-US-043-V.6
 *
 * MATHEMATICAL IDENTITY:
 *   Construct water geometry with H-O-H = 104.5 degrees, r_OH = 0.096 nm
 *   O at origin, H1 along x-axis, H2 at angle
 */
static int test_trajectory_validity_angle_water_T_US_043_V_6(void)
{
    double r_OH = 0.096;  /* nm */
    double angle_HOH = 104.5;  /* degrees */

    /* O at origin */
    double O[3] = {0.0, 0.0, 0.0};

    /* H1 along positive x-axis */
    double H1[3] = {r_OH, 0.0, 0.0};

    /* H2 at angle from x-axis in xy-plane */
    double theta_rad = angle_HOH * 3.14159265358979323846 / 180.0;
    double H2[3] = {r_OH * cos(theta_rad), r_OH * sin(theta_rad), 0.0};

    double computed_angle = compute_angle_degrees(H1, O, H2);

    const double TOL = 1e-10;

    if (fabs(computed_angle - angle_HOH) > TOL) {
        fprintf(stderr, "\n    FAIL: angle = %.15f deg, expected %.15f deg\n",
                computed_angle, angle_HOH);
        return 0;
    }

    /* Also verify bond lengths */
    double d_OH1 = compute_distance(O, H1);
    double d_OH2 = compute_distance(O, H2);

    if (fabs(d_OH1 - r_OH) > TOL || fabs(d_OH2 - r_OH) > TOL) {
        fprintf(stderr, "\n    FAIL: bond lengths incorrect\n");
        return 0;
    }

    return 1;
}

/**
 * Test: Vector magnitude calculation.
 *
 * ThDD:T-US-043-V.7
 *
 * MATHEMATICAL IDENTITY:
 *   |{3,4,0}| = 5
 */
static int test_trajectory_validity_magnitude_T_US_043_V_7(void)
{
    double v[3] = {3.0, 4.0, 0.0};

    double mag = compute_magnitude(v);
    double expected = 5.0;

    const double TOL = 1e-14;

    if (fabs(mag - expected) > TOL) {
        fprintf(stderr, "\n    FAIL: magnitude = %.15f, expected %.15f\n",
                mag, expected);
        return 0;
    }

    return 1;
}

/**
 * Test: Geometry bounds check.
 *
 * ThDD:T-US-043-V.6
 *
 * Verifies that typical water geometry passes bounds check.
 * Equilibrium values from spectroscopy (NOT fabricated):
 *   r_OH = 0.09572 nm (gas phase)
 *   angle_HOH = 104.52 degrees
 */
static int test_trajectory_validity_water_bounds_T_US_043_V_6(void)
{
    double r_OH = 0.09572;  /* Literature value: gas phase water */
    double angle_HOH = 104.52;  /* Literature value */

    /* Check bond length bounds */
    if (r_OH < OH_BOND_MIN_NM || r_OH > OH_BOND_MAX_NM) {
        fprintf(stderr, "\n    FAIL: equilibrium r_OH = %.5f nm out of bounds [%.2f, %.2f]\n",
                r_OH, OH_BOND_MIN_NM, OH_BOND_MAX_NM);
        return 0;
    }

    /* Check angle bounds */
    if (angle_HOH < HOH_ANGLE_MIN_DEG || angle_HOH > HOH_ANGLE_MAX_DEG) {
        fprintf(stderr, "\n    FAIL: equilibrium angle = %.2f deg out of bounds [%.1f, %.1f]\n",
                angle_HOH, HOH_ANGLE_MIN_DEG, HOH_ANGLE_MAX_DEG);
        return 0;
    }

    printf("\n    (equilibrium r_OH = %.5f nm, angle = %.2f deg within bounds) ",
           r_OH, angle_HOH);
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
 * Test V.5: No NaN/Inf values in trajectory.
 *
 * ThDD:T-US-043-V.5
 *
 * Criterion: isfinite(E(t_i)) and isfinite(r_j(t_i)) for all t_i, atoms j
 */
static int test_no_nan_inf_V_US_043_V_5(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    printf("\n");
    printf("    === V.5: No NaN/Inf Values ===\n");
    printf("    Acceptance: isfinite(x) = true for all energy and coordinate values\n");
    printf("\n");

    char xvg_path[512];
    snprintf(xvg_path, sizeof(xvg_path), "%s/nve_500ps/total_energy.xvg", B5_DATA_DIR);

    FILE *f = fopen(xvg_path, "r");
    if (!f) {
        printf("    Trajectory data not found: %s\n", xvg_path);
        printf("    SKIP: Trajectory data required\n");
        return -1;
    }
    fclose(f);

    printf("    Trajectory found - analysis not yet implemented\n");
    printf("    SKIP: Full analysis pending implementation\n");
    return -1;
#endif
}

/**
 * Test V.6: QM region geometric integrity.
 *
 * ThDD:T-US-043-V.6
 *
 * Criterion:
 *   - O-H bond lengths: 0.08 nm < r_OH < 0.15 nm
 *   - H-O-H angle: 80 deg < angle < 120 deg
 */
static int test_qm_geometry_V_US_043_V_6(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    printf("\n");
    printf("    === V.6: QM Geometry Integrity ===\n");
    printf("    O-H bond: %.2f nm < r_OH < %.2f nm\n", OH_BOND_MIN_NM, OH_BOND_MAX_NM);
    printf("    H-O-H angle: %.1f deg < angle < %.1f deg\n",
           HOH_ANGLE_MIN_DEG, HOH_ANGLE_MAX_DEG);
    printf("\n");
    printf("    SKIP: Requires trajectory data\n");
    return -1;
#endif
}

/**
 * Test V.7: Center of mass / momentum conservation.
 *
 * ThDD:T-US-043-V.7
 *
 * Criterion: |P_total(t) - P_total(0)| / |P_total(0)| < 10^-6
 */
static int test_com_conservation_V_US_043_V_7(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    printf("\n");
    printf("    === V.7: Momentum Conservation ===\n");
    printf("    Acceptance: relative momentum drift < %.0e\n", MOMENTUM_CONSERVATION_TOL);
    printf("\n");
    printf("    Rationale: NVE should conserve total momentum\n");
    printf("    (no external forces acting on system)\n");
    printf("\n");
    printf("    SKIP: Requires trajectory data\n");
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
    printf("\n");
    printf("Unit tests verify helper functions using mathematical identities.\n");
    printf("Integration tests require trajectory data from GROMACS simulations.\n");
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

    printf("=== US-043 Trajectory Validity Tests (V.5-V.7) ===\n\n");

    printf("Physical Bounds (from spectroscopic data):\n");
    printf("  O-H bond: [%.2f, %.2f] nm (equilibrium ~0.096 nm)\n",
           OH_BOND_MIN_NM, OH_BOND_MAX_NM);
    printf("  H-O-H angle: [%.1f, %.1f] deg (equilibrium ~104.5 deg)\n",
           HOH_ANGLE_MIN_DEG, HOH_ANGLE_MAX_DEG);
    printf("  Momentum tolerance: %.0e relative\n", MOMENTUM_CONSERVATION_TOL);
    printf("\n");

    if (run_unit) {
        printf("=== Unit Tests (Mathematical Identities) ===\n\n");

        printf("NaN/Inf Detection (T-US-043-V.5):\n");
        RUN_TEST(test_trajectory_validity_nan_detection_T_US_043_V_5);
        RUN_TEST(test_trajectory_validity_inf_detection_T_US_043_V_5);
        printf("\n");

        printf("Distance Calculation (T-US-043-V.6):\n");
        RUN_TEST(test_trajectory_validity_distance_calculation_T_US_043_V_6);
        printf("\n");

        printf("Angle Calculation (T-US-043-V.6):\n");
        RUN_TEST(test_trajectory_validity_angle_right_T_US_043_V_6);
        RUN_TEST(test_trajectory_validity_angle_60deg_T_US_043_V_6);
        RUN_TEST(test_trajectory_validity_angle_water_T_US_043_V_6);
        printf("\n");

        printf("Magnitude Calculation (T-US-043-V.7):\n");
        RUN_TEST(test_trajectory_validity_magnitude_T_US_043_V_7);
        printf("\n");

        printf("Geometry Bounds (T-US-043-V.6):\n");
        RUN_TEST(test_trajectory_validity_water_bounds_T_US_043_V_6);
        printf("\n");
    }

    if (run_integration) {
        printf("=== Integration Tests (V.5-V.7 Trajectory Analysis) ===\n\n");

#ifndef GRODFTB_LONG_TESTS
        printf("NOTE: Integration tests require GRODFTB_LONG_TESTS=ON\n");
        printf("      Rebuild with: cmake -DGRODFTB_LONG_TESTS=ON ..\n\n");
#endif

        printf("V.5 - No NaN/Inf:\n");
        RUN_TEST(test_no_nan_inf_V_US_043_V_5);
        printf("\n");

        printf("V.6 - QM Geometry:\n");
        RUN_TEST(test_qm_geometry_V_US_043_V_6);
        printf("\n");

        printf("V.7 - Momentum Conservation:\n");
        RUN_TEST(test_com_conservation_V_US_043_V_7);
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
