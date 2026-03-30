/*
 * US-043: Checkpoint/Restart Verification Tests (V.13-V.14)
 *
 * ThDD:T-US-043-V.13 -- Checkpoint restart determinism
 * ThDD:T-US-043-V.14 -- Energy continuity at restart
 * SDD:specs.md:S14 -- Checkpoint and Restart interface specification
 * SDD:docs/verification/US-043.md -- Verification plan
 *
 * This file implements:
 *   1. V.13: Restart determinism - trajectory restarted from checkpoint must
 *            produce identical results to continuous trajectory
 *   2. V.14: Energy continuity - no energy jump at restart point
 *
 * GOLDEN RULE COMPLIANCE:
 * - Unit tests use mathematical identities (exact formulas)
 * - Acceptance thresholds from specification (10^-10 relative tolerance)
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
#include <stdint.h>

#include "grodftb/units.h"
#include "grodftb/error.h"

/* ===========================================================================
 * ACCEPTANCE CRITERIA (from T-US-043-V.13 and V.14)
 *
 * SDD:specs.md:S14.1 -- Checkpoint restart determinism requirement:
 *   "...so that a restarted simulation produces identical results to an
 *    uninterrupted one."
 *
 * ThDD:T-US-043-V.13:
 *   |E_restart(t) - E_continuous(t)| < 10^-10 * |E|
 *
 * ThDD:T-US-043-V.14:
 *   |E(t_checkpoint+) - E(t_checkpoint-)| < 10^-10 * |E|
 *
 * These are specification-derived thresholds, NOT fabricated expected values.
 * The 10^-10 relative tolerance is chosen to allow for:
 *   - Double precision round-off (~10^-16)
 *   - Accumulated error over many operations (~10^-12 to 10^-10)
 *   - But NOT algorithmic differences or missing state
 * ===========================================================================*/

/* ThDD:T-US-043-V.13, V.14 -- Relative tolerance for restart determinism */
#define CHECKPOINT_RELATIVE_TOL  1e-10

/* Absolute tolerance for energies near zero (to avoid division by zero) */
#define CHECKPOINT_ABSOLUTE_TOL  1e-20

/* Maximum number of steps to compare after restart */
#define MAX_COMPARE_STEPS  10000

/* ===========================================================================
 * B5 SYSTEM PARAMETERS (from tests/data/b5/provenance.json)
 *
 * ThDD:CLAUDE.md:Golden_Rule -- All values from actual GROMACS preparation
 * ===========================================================================*/

#define B5_N_ATOMS_TOTAL  2652   /* Total atoms in system */
#define B5_N_QM_ATOMS     3      /* QM region: 1 water molecule */

/* ===========================================================================
 * SECTION 1: MATHEMATICAL HELPER FUNCTIONS
 *
 * ThDD: These implement the relative error calculation specified in V.13/V.14.
 * All expected values are mathematical identities.
 * ===========================================================================*/

/**
 * Compute relative error between two values.
 *
 * ThDD:T-US-043-V.13 -- Relative error formula
 *
 * relative_error = |a - b| / max(|a|, |b|, eps)
 *
 * The denominator uses the maximum of the two absolute values to provide
 * a symmetric measure. A small epsilon prevents division by zero when
 * both values are near zero.
 *
 * @param a       First value
 * @param b       Second value
 * @param eps     Minimum denominator (prevents division by zero)
 * @return        Relative error (dimensionless)
 */
static double compute_relative_error(double a, double b, double eps)
{
    double diff = fabs(a - b);
    double scale = fabs(a);
    if (fabs(b) > scale) {
        scale = fabs(b);
    }
    if (scale < eps) {
        scale = eps;
    }
    return diff / scale;
}

/**
 * Compute absolute error between two values.
 *
 * ThDD:T-US-043-V.14 -- Absolute difference
 *
 * @param a  First value
 * @param b  Second value
 * @return   Absolute difference |a - b|
 */
static double __attribute__((unused)) compute_absolute_error(double a, double b)
{
    return fabs(a - b);
}

/**
 * Check if two energy values match within checkpoint tolerance.
 *
 * ThDD:T-US-043-V.13, V.14 -- Combined relative and absolute check
 *
 * Uses BOTH relative and absolute tolerances:
 *   PASS if |a - b| < max(rel_tol * max(|a|,|b|), abs_tol)
 *
 * This handles both large energies (where relative error matters) and
 * small energies near zero (where absolute error matters).
 *
 * @param a        First energy value
 * @param b        Second energy value
 * @param rel_tol  Relative tolerance (e.g., 1e-10)
 * @param abs_tol  Absolute tolerance (e.g., 1e-20)
 * @return         1 if values match within tolerance, 0 otherwise
 */
static int energies_match(double a, double b, double rel_tol, double abs_tol)
{
    double diff = fabs(a - b);
    double scale = fabs(a);
    if (fabs(b) > scale) {
        scale = fabs(b);
    }

    double threshold = rel_tol * scale;
    if (threshold < abs_tol) {
        threshold = abs_tol;
    }

    return (diff <= threshold) ? 1 : 0;
}

/**
 * Find the maximum relative error in a trajectory comparison.
 *
 * ThDD:T-US-043-V.13 -- Max relative error over trajectory segment
 *
 * @param n           Number of points
 * @param E_cont      Continuous trajectory energies [n]
 * @param E_restart   Restarted trajectory energies [n]
 * @param eps         Minimum scale for relative error
 * @param max_idx_out Receives index of maximum error (may be NULL)
 * @return            Maximum relative error
 */
static double find_max_relative_error(int n, const double *E_cont,
                                       const double *E_restart, double eps,
                                       int *max_idx_out)
{
    double max_err = 0.0;
    int max_idx = 0;

    for (int i = 0; i < n; i++) {
        double rel_err = compute_relative_error(E_cont[i], E_restart[i], eps);
        if (rel_err > max_err) {
            max_err = rel_err;
            max_idx = i;
        }
    }

    if (max_idx_out) {
        *max_idx_out = max_idx;
    }

    return max_err;
}

/**
 * Compute RMS relative error over a trajectory segment.
 *
 * ThDD:T-US-043-V.13 -- RMS error diagnostic
 *
 * @param n          Number of points
 * @param E_cont     Continuous trajectory energies [n]
 * @param E_restart  Restarted trajectory energies [n]
 * @param eps        Minimum scale for relative error
 * @return           RMS relative error
 */
static double compute_rms_relative_error(int n, const double *E_cont,
                                          const double *E_restart, double eps)
{
    if (n <= 0) return 0.0;

    double sum_sq = 0.0;
    for (int i = 0; i < n; i++) {
        double rel_err = compute_relative_error(E_cont[i], E_restart[i], eps);
        sum_sq += rel_err * rel_err;
    }

    return sqrt(sum_sq / (double)n);
}

/* ===========================================================================
 * SECTION 2: CHECKPOINT BUFFER VALIDATION HELPERS
 *
 * ThDD:T-US-043-V.13 -- State completeness verification
 * SDD:specs.md:S14.2 -- State that must be serialized
 *
 * These functions help verify checkpoint buffer integrity without
 * fabricating expected values.
 * ===========================================================================*/

/**
 * Checkpoint header structure for validation.
 *
 * SDD:specs.md:S14.2 -- Minimum state components
 *
 * This is a conceptual structure for testing the checkpoint format.
 * The actual implementation may differ.
 */
typedef struct {
    uint32_t magic;           /* Magic number for format identification */
    uint32_t version;         /* Checkpoint format version */
    uint32_t n_qm_atoms;      /* Number of QM atoms */
    uint32_t has_scc_charges; /* Flag: SCC charges present */
    uint32_t has_fssh_state;  /* Flag: FSSH state present */
    uint32_t has_ehrenfest;   /* Flag: Ehrenfest state present */
    uint32_t data_size;       /* Size of payload data in bytes */
} checkpoint_header_t;

/* Expected magic number for GRO-DFTB checkpoint format */
#define GRODFTB_CHECKPOINT_MAGIC  0x47524F44  /* "GROD" in ASCII */

/**
 * Validate checkpoint buffer has minimum required structure.
 *
 * SDD:specs.md:S14.3 -- Serialization interface
 *
 * This does NOT validate content (which would require fabricating expected
 * values), only that the buffer has plausible structure.
 *
 * @param buf   Checkpoint buffer
 * @param size  Buffer size in bytes
 * @return      1 if structure appears valid, 0 otherwise
 */
static int __attribute__((unused)) validate_checkpoint_structure(const void *buf, size_t size)
{
    if (!buf || size < sizeof(checkpoint_header_t)) {
        return 0;
    }

    const checkpoint_header_t *hdr = (const checkpoint_header_t *)buf;

    /* Verify magic number */
    if (hdr->magic != GRODFTB_CHECKPOINT_MAGIC) {
        return 0;
    }

    /* Verify size consistency */
    if (size < sizeof(checkpoint_header_t) + hdr->data_size) {
        return 0;
    }

    return 1;
}

/**
 * Compute checksum of checkpoint buffer (for identity testing).
 *
 * ThDD:T-US-043-V.13 -- Deterministic state verification
 *
 * Uses a simple FNV-1a hash for detecting differences.
 *
 * @param buf   Checkpoint buffer
 * @param size  Buffer size in bytes
 * @return      Hash value
 */
static uint64_t compute_checkpoint_hash(const void *buf, size_t size)
{
    /*
     * FNV-1a hash for 64-bit
     * ThDD: Mathematical identity - same input always produces same hash
     */
    const uint64_t FNV_OFFSET = 14695981039346656037ULL;
    const uint64_t FNV_PRIME = 1099511628211ULL;

    uint64_t hash = FNV_OFFSET;
    const uint8_t *bytes = (const uint8_t *)buf;

    for (size_t i = 0; i < size; i++) {
        hash ^= bytes[i];
        hash *= FNV_PRIME;
    }

    return hash;
}

/* ===========================================================================
 * SECTION 3: TEST FRAMEWORK
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
 * SECTION 4: UNIT TESTS (Mathematical Identities)
 *
 * GOLDEN RULE: All expected values are mathematical identities,
 * not fabricated physical results.
 * ===========================================================================*/

/**
 * Test: Relative error calculation for identical values.
 *
 * MATHEMATICAL IDENTITY:
 *   relative_error(x, x) = 0 for any x
 */
static int test_checkpoint_relative_error_identical_T_US_043_N_1(void)
{
    double values[] = {1.0, 100.0, -50.0, 1e10, 1e-10, 0.0};
    int n = sizeof(values) / sizeof(values[0]);

    for (int i = 0; i < n; i++) {
        double err = compute_relative_error(values[i], values[i],
                                            CHECKPOINT_ABSOLUTE_TOL);
        if (err != 0.0) {
            fprintf(stderr, "\n    FAIL: relative_error(%g, %g) = %g, expected 0\n",
                    values[i], values[i], err);
            return 0;
        }
    }

    return 1;
}

/**
 * Test: Relative error calculation for known differences.
 *
 * MATHEMATICAL IDENTITY:
 *   relative_error(100, 101) = |100-101| / max(100, 101) = 1/101 ~ 0.0099
 *   relative_error(100, 99)  = |100-99| / max(100, 99)   = 1/100 = 0.01
 */
static int test_checkpoint_relative_error_known_T_US_043_N_2(void)
{
    const double TOL = 1e-14;

    /* Test case 1: 100 vs 101 */
    double err1 = compute_relative_error(100.0, 101.0, CHECKPOINT_ABSOLUTE_TOL);
    double expected1 = 1.0 / 101.0;

    if (fabs(err1 - expected1) > TOL) {
        fprintf(stderr, "\n    FAIL: relative_error(100, 101) = %.15e, expected %.15e\n",
                err1, expected1);
        return 0;
    }

    /* Test case 2: 100 vs 99 */
    double err2 = compute_relative_error(100.0, 99.0, CHECKPOINT_ABSOLUTE_TOL);
    double expected2 = 1.0 / 100.0;

    if (fabs(err2 - expected2) > TOL) {
        fprintf(stderr, "\n    FAIL: relative_error(100, 99) = %.15e, expected %.15e\n",
                err2, expected2);
        return 0;
    }

    /* Test case 3: Symmetry check - error(a,b) == error(b,a) */
    double err3a = compute_relative_error(100.0, 105.0, CHECKPOINT_ABSOLUTE_TOL);
    double err3b = compute_relative_error(105.0, 100.0, CHECKPOINT_ABSOLUTE_TOL);

    if (fabs(err3a - err3b) > TOL) {
        fprintf(stderr, "\n    FAIL: Asymmetric error: err(100,105)=%g, err(105,100)=%g\n",
                err3a, err3b);
        return 0;
    }

    return 1;
}

/**
 * Test: Relative error with values near zero.
 *
 * ThDD:T-US-043-V.14 -- Handling small energies
 *
 * MATHEMATICAL IDENTITY:
 *   When |a| and |b| are both < eps, the scale becomes eps:
 *   relative_error(1e-30, 2e-30) with eps=1e-20 = |1e-30 - 2e-30| / 1e-20 = 1e-10
 */
static int test_checkpoint_relative_error_near_zero_T_US_043_N_3(void)
{
    const double eps = 1e-20;
    const double TOL = 1e-14;

    double a = 1e-30;
    double b = 2e-30;

    double err = compute_relative_error(a, b, eps);
    double expected = fabs(a - b) / eps;  /* = 1e-30 / 1e-20 = 1e-10 */

    if (fabs(err - expected) > TOL) {
        fprintf(stderr, "\n    FAIL: relative_error(1e-30, 2e-30) = %.15e, expected %.15e\n",
                err, expected);
        return 0;
    }

    printf("\n    (small values handled correctly: err = %.2e) ", err);
    return 1;
}

/**
 * Test: energies_match function for identical values.
 *
 * MATHEMATICAL IDENTITY:
 *   energies_match(x, x) = true for any x (within any positive tolerance)
 */
static int test_checkpoint_energies_match_identical_T_US_043_N_4(void)
{
    double values[] = {-1000.0, 0.0, 1e-20, 1e20};
    int n = sizeof(values) / sizeof(values[0]);

    for (int i = 0; i < n; i++) {
        if (!energies_match(values[i], values[i],
                            CHECKPOINT_RELATIVE_TOL, CHECKPOINT_ABSOLUTE_TOL)) {
            fprintf(stderr, "\n    FAIL: energies_match(%g, %g) returned false\n",
                    values[i], values[i]);
            return 0;
        }
    }

    return 1;
}

/**
 * Test: energies_match function for values within tolerance.
 *
 * ThDD:T-US-043-V.13, V.14 -- Tolerance boundary behavior
 *
 * MATHEMATICAL IDENTITY:
 *   With rel_tol = 1e-10, energies_match(1e6, 1e6*(1 + 0.5e-10)) should be true
 *   because the relative error is 0.5e-10 < 1e-10
 */
static int test_checkpoint_energies_match_within_tol_T_US_043_N_5(void)
{
    double base = 1e6;  /* Representative energy magnitude (kJ/mol) */

    /* Value within 50% of tolerance should match */
    double diff_ok = base * (0.5 * CHECKPOINT_RELATIVE_TOL);
    if (!energies_match(base, base + diff_ok,
                        CHECKPOINT_RELATIVE_TOL, CHECKPOINT_ABSOLUTE_TOL)) {
        fprintf(stderr, "\n    FAIL: Values within 50%% tolerance should match\n");
        fprintf(stderr, "    base = %g, base + diff = %g\n", base, base + diff_ok);
        return 0;
    }

    /* Value at 200% of tolerance should NOT match */
    double diff_bad = base * (2.0 * CHECKPOINT_RELATIVE_TOL);
    if (energies_match(base, base + diff_bad,
                       CHECKPOINT_RELATIVE_TOL, CHECKPOINT_ABSOLUTE_TOL)) {
        fprintf(stderr, "\n    FAIL: Values at 200%% tolerance should NOT match\n");
        fprintf(stderr, "    base = %g, base + diff = %g\n", base, base + diff_bad);
        return 0;
    }

    printf("\n    (boundary behavior verified for tol = %.0e) ", CHECKPOINT_RELATIVE_TOL);
    return 1;
}

/**
 * Test: Maximum relative error finder.
 *
 * MATHEMATICAL IDENTITY:
 *   For arrays differing only at index 5 by factor (1 + 1e-8):
 *   max_relative_error should be ~1e-8 and occur at index 5
 */
static int test_checkpoint_max_relative_error_T_US_043_N_6(void)
{
    const int n = 10;
    double E_cont[10] = {-100.0, -100.1, -100.2, -100.3, -100.4,
                         -100.5, -100.6, -100.7, -100.8, -100.9};
    double E_restart[10];

    /* Copy and introduce a single difference at index 5 */
    for (int i = 0; i < n; i++) {
        E_restart[i] = E_cont[i];
    }
    E_restart[5] = E_cont[5] * (1.0 + 1e-8);

    int max_idx = -1;
    double max_err = find_max_relative_error(n, E_cont, E_restart,
                                              CHECKPOINT_ABSOLUTE_TOL, &max_idx);

    /* Should find the error at index 5 */
    if (max_idx != 5) {
        fprintf(stderr, "\n    FAIL: max error at index %d, expected 5\n", max_idx);
        return 0;
    }

    /* Error magnitude should be ~1e-8 */
    if (fabs(max_err - 1e-8) > 1e-10) {
        fprintf(stderr, "\n    FAIL: max_err = %.2e, expected ~1e-8\n", max_err);
        return 0;
    }

    return 1;
}

/**
 * Test: RMS relative error for identical arrays.
 *
 * MATHEMATICAL IDENTITY:
 *   RMS error of (x, x) arrays = 0
 */
static int test_checkpoint_rms_error_identical_T_US_043_N_7(void)
{
    const int n = 100;
    double *E = malloc(n * sizeof(double));
    if (!E) {
        fprintf(stderr, "\n    FAIL: Memory allocation failed\n");
        return 0;
    }

    for (int i = 0; i < n; i++) {
        E[i] = -1000.0 + i * 0.1;
    }

    double rms = compute_rms_relative_error(n, E, E, CHECKPOINT_ABSOLUTE_TOL);
    free(E);

    if (rms != 0.0) {
        fprintf(stderr, "\n    FAIL: RMS of identical arrays = %g, expected 0\n", rms);
        return 0;
    }

    return 1;
}

/**
 * Test: FNV-1a hash determinism.
 *
 * MATHEMATICAL IDENTITY:
 *   hash(x) == hash(x) for same input
 *   hash(x) != hash(y) for different inputs (with high probability)
 */
static int test_checkpoint_hash_determinism_T_US_043_N_8(void)
{
    uint8_t data1[] = {0x01, 0x02, 0x03, 0x04, 0x05};
    uint8_t data2[] = {0x01, 0x02, 0x03, 0x04, 0x06};  /* Last byte different */

    uint64_t hash1a = compute_checkpoint_hash(data1, sizeof(data1));
    uint64_t hash1b = compute_checkpoint_hash(data1, sizeof(data1));
    uint64_t hash2 = compute_checkpoint_hash(data2, sizeof(data2));

    /* Same data should give same hash */
    if (hash1a != hash1b) {
        fprintf(stderr, "\n    FAIL: Non-deterministic hash\n");
        return 0;
    }

    /* Different data should give different hash (extremely high probability) */
    if (hash1a == hash2) {
        fprintf(stderr, "\n    FAIL: Hash collision on different data\n");
        return 0;
    }

    return 1;
}

/**
 * Test: Checkpoint tolerance is within machine precision limits.
 *
 * ThDD:T-US-043-V.13 -- Tolerance justification
 *
 * MATHEMATICAL IDENTITY:
 *   The specified tolerance (1e-10) must be achievable with double precision:
 *   - Machine epsilon for double: ~2.2e-16
 *   - Practical accumulated error: ~1e-12 to 1e-10
 *   - Therefore 1e-10 is achievable but stringent
 */
static int test_checkpoint_tolerance_achievable_T_US_043_N_9(void)
{
    /* Verify tolerance is stricter than single precision epsilon (~1e-7)
     * but looser than double precision epsilon (~2e-16) */
    if (CHECKPOINT_RELATIVE_TOL >= 1e-7) {
        fprintf(stderr, "\n    FAIL: Tolerance %.0e is too loose (single precision)\n",
                CHECKPOINT_RELATIVE_TOL);
        return 0;
    }

    if (CHECKPOINT_RELATIVE_TOL < DBL_EPSILON * 1e4) {
        fprintf(stderr, "\n    FAIL: Tolerance %.0e may be unachievable\n",
                CHECKPOINT_RELATIVE_TOL);
        fprintf(stderr, "    (< 10000 * DBL_EPSILON = %.2e)\n", DBL_EPSILON * 1e4);
        return 0;
    }

    printf("\n    (tol = %.0e, DBL_EPSILON = %.2e) ", CHECKPOINT_RELATIVE_TOL, DBL_EPSILON);
    return 1;
}

/* ===========================================================================
 * SECTION 5: INTEGRATION TESTS (Trajectory Comparison)
 *
 * These tests require actual trajectory data from GROMACS simulations:
 *   1. A continuous trajectory (no restart)
 *   2. A trajectory restarted from a checkpoint at time t_checkpoint
 *
 * Gated behind GRODFTB_LONG_TESTS.
 *
 * GOLDEN RULE: All expected values are thresholds from specs.md, NOT
 * fabricated energy predictions. We compare trajectory to trajectory,
 * not trajectory to pre-computed values.
 * ===========================================================================*/

#ifndef GRODFTB_SOURCE_DIR_STR
#define GRODFTB_SOURCE_DIR_STR "."
#endif
#define B5_DATA_DIR GRODFTB_SOURCE_DIR_STR "/tests/data/b5"

/**
 * Test V.13: Checkpoint restart determinism.
 *
 * ThDD:T-US-043-V.13
 * SDD:specs.md:S14.1
 *
 * Criterion: |E_restart(t) - E_continuous(t)| < 10^-10 * |E|
 *            for all t > t_checkpoint
 *
 * METHODOLOGY:
 * 1. Run continuous trajectory: 0 -> T (e.g., 10 ps)
 * 2. Run trajectory with restart: 0 -> T/2 (checkpoint), restart -> T
 * 3. Compare energies after restart point
 * 4. Maximum relative error must be < 10^-10
 *
 * NOTE: This test SKIPS if trajectory data is not available.
 * The trajectories must be generated by running GROMACS mdrun with
 * appropriate checkpoint settings.
 */
static int test_restart_determinism_V_US_043_V_13(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    printf("\n");
    printf("    === V.13: Checkpoint Restart Determinism ===\n");
    printf("    Criterion: max |E_restart - E_cont| / |E| < %.0e\n",
           CHECKPOINT_RELATIVE_TOL);
    printf("    Source: specs.md S14.1, T-US-043-V.13\n");
    printf("\n");

    /* Check for continuous trajectory */
    char cont_path[512];
    snprintf(cont_path, sizeof(cont_path),
             "%s/checkpoint_test/continuous_energy.xvg", B5_DATA_DIR);

    FILE *f_cont = fopen(cont_path, "r");
    if (!f_cont) {
        printf("    Continuous trajectory not found:\n");
        printf("      %s\n", cont_path);
        printf("\n");
        printf("    To generate test data:\n");
        printf("    1. Run continuous simulation:\n");
        printf("       gmx mdrun -s nve.tpr -deffnm continuous\n");
        printf("    2. Run simulation with checkpoint at 5ps:\n");
        printf("       gmx mdrun -s nve.tpr -deffnm checkpoint -cpt 5\n");
        printf("    3. Kill and restart from checkpoint:\n");
        printf("       gmx mdrun -s nve.tpr -deffnm restart -cpi checkpoint.cpt\n");
        printf("    4. Extract energies:\n");
        printf("       gmx energy -f continuous.edr -o continuous_energy.xvg\n");
        printf("       gmx energy -f restart.edr -o restart_energy.xvg\n");
        printf("\n");
        printf("    SKIP: Trajectory data required\n");
        return -1;
    }
    fclose(f_cont);

    /* Check for restarted trajectory */
    char restart_path[512];
    snprintf(restart_path, sizeof(restart_path),
             "%s/checkpoint_test/restart_energy.xvg", B5_DATA_DIR);

    FILE *f_restart = fopen(restart_path, "r");
    if (!f_restart) {
        printf("    Restarted trajectory not found:\n");
        printf("      %s\n", restart_path);
        printf("\n    SKIP: Trajectory data required\n");
        return -1;
    }
    fclose(f_restart);

    /* TODO: Implement XVG parsing and trajectory comparison */
    printf("    Both trajectory files found.\n");
    printf("    Analysis not yet implemented.\n");
    printf("\n");
    printf("    Expected analysis:\n");
    printf("    - Parse time-energy pairs from both files\n");
    printf("    - Align time points after checkpoint\n");
    printf("    - Compute relative error at each point\n");
    printf("    - Report max error and PASS/FAIL\n");
    printf("\n    SKIP: Full analysis pending implementation\n");
    return -1;
#endif
}

/**
 * Test V.14: Energy continuity at restart point.
 *
 * ThDD:T-US-043-V.14
 * SDD:specs.md:S14.1
 *
 * Criterion: |E(t_checkpoint+) - E(t_checkpoint-)| < 10^-10 * |E|
 *
 * METHODOLOGY:
 * 1. Extract energy immediately before checkpoint (from continuous run)
 * 2. Extract energy immediately after restart (from restarted run)
 * 3. These should be essentially identical (within numerical precision)
 *
 * A jump in energy at restart indicates incomplete state restoration.
 *
 * NOTE: This test SKIPS if trajectory data is not available.
 */
static int test_restart_energy_continuity_V_US_043_V_14(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    printf("\n");
    printf("    === V.14: Energy Continuity at Restart ===\n");
    printf("    Criterion: |E(t+) - E(t-)| / |E| < %.0e at checkpoint\n",
           CHECKPOINT_RELATIVE_TOL);
    printf("    Source: specs.md S14.1, T-US-043-V.14\n");
    printf("\n");

    /* Check for checkpoint-time energy files */
    char before_path[512];
    snprintf(before_path, sizeof(before_path),
             "%s/checkpoint_test/energy_at_checkpoint.dat", B5_DATA_DIR);

    FILE *f_before = fopen(before_path, "r");
    if (!f_before) {
        printf("    Energy data at checkpoint not found:\n");
        printf("      %s\n", before_path);
        printf("\n");
        printf("    To generate test data:\n");
        printf("    1. Identify checkpoint time (e.g., t_cpt = 5.0 ps)\n");
        printf("    2. Extract energy at t_cpt from continuous run:\n");
        printf("       gmx energy -f continuous.edr -o energy_at_checkpoint.dat -b 5.0 -e 5.001\n");
        printf("    3. Extract first energy after restart:\n");
        printf("       gmx energy -f restart.edr -o energy_after_restart.dat -b 5.0 -e 5.001\n");
        printf("\n");
        printf("    SKIP: Checkpoint energy data required\n");
        return -1;
    }
    fclose(f_before);

    char after_path[512];
    snprintf(after_path, sizeof(after_path),
             "%s/checkpoint_test/energy_after_restart.dat", B5_DATA_DIR);

    FILE *f_after = fopen(after_path, "r");
    if (!f_after) {
        printf("    Energy after restart not found:\n");
        printf("      %s\n", after_path);
        printf("\n    SKIP: Checkpoint energy data required\n");
        return -1;
    }
    fclose(f_after);

    /* TODO: Implement energy file parsing and comparison */
    printf("    Both energy files found.\n");
    printf("    Analysis not yet implemented.\n");
    printf("\n");
    printf("    Expected analysis:\n");
    printf("    - Read energy at t_checkpoint from continuous run\n");
    printf("    - Read first energy from restarted run (at t_checkpoint+)\n");
    printf("    - Compute relative difference\n");
    printf("    - Compare against tolerance (%.0e)\n", CHECKPOINT_RELATIVE_TOL);
    printf("    - An energy jump indicates incomplete checkpoint state\n");
    printf("\n    SKIP: Full analysis pending implementation\n");
    return -1;
#endif
}

/**
 * Test: Multi-restart determinism.
 *
 * ThDD:T-US-043-V.13 (extended)
 *
 * Verifies that multiple sequential restarts don't accumulate error.
 *
 * METHODOLOGY:
 * 1. Run: 0 -> T/3 (cpt1) -> restart -> 2T/3 (cpt2) -> restart -> T
 * 2. Compare against continuous 0 -> T
 * 3. Error after second restart should still be < 10^-10
 *
 * NOTE: This test SKIPS if trajectory data is not available.
 */
static int test_multi_restart_determinism_V_US_043_V_13_ext(void)
{
#ifndef GRODFTB_LONG_TESTS
    printf("\n    SKIP: Long test (requires GRODFTB_LONG_TESTS)\n");
    return -1;
#else
    printf("\n");
    printf("    === V.13 Extended: Multi-Restart Determinism ===\n");
    printf("    Criterion: Error after multiple restarts still < %.0e\n",
           CHECKPOINT_RELATIVE_TOL);
    printf("\n");

    char multi_path[512];
    snprintf(multi_path, sizeof(multi_path),
             "%s/checkpoint_test/multi_restart_energy.xvg", B5_DATA_DIR);

    FILE *f = fopen(multi_path, "r");
    if (!f) {
        printf("    Multi-restart trajectory not found:\n");
        printf("      %s\n", multi_path);
        printf("\n");
        printf("    This is an extended validation test.\n");
        printf("    Generate by running simulation with multiple restarts.\n");
        printf("\n    SKIP: Multi-restart data not available\n");
        return -1;
    }
    fclose(f);

    printf("    Data found - analysis not yet implemented\n");
    printf("    SKIP: Full analysis pending implementation\n");
    return -1;
#endif
}

/* ===========================================================================
 * SECTION 6: MAIN TEST RUNNER
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
    printf("Checkpoint/Restart Verification (V.13-V.14):\n");
    printf("  V.13: Restart determinism (trajectory comparison)\n");
    printf("  V.14: Energy continuity at restart point\n");
    printf("\n");
    printf("Unit tests verify mathematical helper functions.\n");
    printf("Integration tests require trajectory data from paired GROMACS runs:\n");
    printf("  1. Continuous trajectory (no restart)\n");
    printf("  2. Trajectory restarted from checkpoint\n");
    printf("\n");
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

    printf("=== US-043 Checkpoint/Restart Verification Tests (V.13-V.14) ===\n\n");

    printf("Acceptance Criteria (from specs.md S14, T-US-043-V.13, V.14):\n");
    printf("  V.13: Restart determinism:    |E_restart - E_cont| / |E| < %.0e\n",
           CHECKPOINT_RELATIVE_TOL);
    printf("  V.14: Energy continuity:      |E(t+) - E(t-)| / |E| < %.0e at checkpoint\n",
           CHECKPOINT_RELATIVE_TOL);
    printf("\n");

    printf("Tolerance Justification (T-US-043-V.13):\n");
    printf("  - Machine precision (double): ~ %.0e\n", DBL_EPSILON);
    printf("  - Accumulated numerical error: ~ 1e-12 to 1e-10\n");
    printf("  - Specified tolerance: %.0e (achievable with exact state restoration)\n",
           CHECKPOINT_RELATIVE_TOL);
    printf("  - Tolerance from specs.md S14.1, verified against numerical analysis\n");
    printf("\n");

    printf("State to Serialize (specs.md S14.2):\n");
    printf("  - SCC converged charges (dq)\n");
    printf("  - Wavefunction coefficients (if excited-state)\n");
    printf("  - FSSH active state and electronic coefficients (if surface hopping)\n");
    printf("  - Ehrenfest coefficients (if Ehrenfest dynamics)\n");
    printf("  - RNG state for stochastic algorithms\n");
    printf("\n");

    if (run_unit) {
        printf("=== Unit Tests (Mathematical Identities) ===\n\n");

        printf("Relative Error Calculation:\n");
        RUN_TEST(test_checkpoint_relative_error_identical_T_US_043_N_1);
        RUN_TEST(test_checkpoint_relative_error_known_T_US_043_N_2);
        RUN_TEST(test_checkpoint_relative_error_near_zero_T_US_043_N_3);
        printf("\n");

        printf("Energy Matching Function:\n");
        RUN_TEST(test_checkpoint_energies_match_identical_T_US_043_N_4);
        RUN_TEST(test_checkpoint_energies_match_within_tol_T_US_043_N_5);
        printf("\n");

        printf("Trajectory Comparison Helpers:\n");
        RUN_TEST(test_checkpoint_max_relative_error_T_US_043_N_6);
        RUN_TEST(test_checkpoint_rms_error_identical_T_US_043_N_7);
        printf("\n");

        printf("Checkpoint Buffer Helpers:\n");
        RUN_TEST(test_checkpoint_hash_determinism_T_US_043_N_8);
        printf("\n");

        printf("Tolerance Validation:\n");
        RUN_TEST(test_checkpoint_tolerance_achievable_T_US_043_N_9);
        printf("\n");
    }

    if (run_integration) {
        printf("=== Integration Tests (V.13-V.14 Trajectory Comparison) ===\n\n");

#ifndef GRODFTB_LONG_TESTS
        printf("NOTE: Integration tests require GRODFTB_LONG_TESTS=ON\n");
        printf("      Rebuild with: cmake -DGRODFTB_LONG_TESTS=ON ..\n\n");
#endif

        printf("V.13 - Restart Determinism:\n");
        RUN_TEST(test_restart_determinism_V_US_043_V_13);
        printf("\n");

        printf("V.14 - Energy Continuity:\n");
        RUN_TEST(test_restart_energy_continuity_V_US_043_V_14);
        printf("\n");

        printf("V.13 Extended - Multi-Restart:\n");
        RUN_TEST(test_multi_restart_determinism_V_US_043_V_13_ext);
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
        printf("  2. Paired trajectory data in tests/data/b5/checkpoint_test/:\n");
        printf("     - continuous_energy.xvg (uninterrupted trajectory)\n");
        printf("     - restart_energy.xvg (trajectory with checkpoint restart)\n");
        printf("     - energy_at_checkpoint.dat (energy just before checkpoint)\n");
        printf("     - energy_after_restart.dat (energy just after restart)\n");
        printf("\n");
        printf("  Generate data by running GROMACS with checkpoint enabled:\n");
        printf("    gmx mdrun -s nve.tpr -cpt 5  # Checkpoint every 5 ps\n");
    }

    if (tests_failed > 0) {
        printf("\nFAILED: %d test(s) failed\n", tests_failed);
        return 1;
    }

    return 0;
}
