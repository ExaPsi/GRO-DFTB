/*
 * TDD:US-006 — Unit tests for CODATA 2022 constants and conversion functions.
 * SDD:specs.md:§6.1–§6.3
 *
 * Tests: 11 acceptance criteria from US-006.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "grodftb/units.h"

static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

#define ASSERT_DOUBLE_EQ(name, a, b)                                         \
    do {                                                                      \
        double _a = (a), _b = (b);                                           \
        if (_a == _b) {                                                       \
            tests_passed++;                                                   \
        } else {                                                              \
            fprintf(stderr, "FAIL %s: %.17g != %.17g\n", name, _a, _b);      \
            tests_failed++;                                                   \
        }                                                                     \
        tests_run++;                                                          \
    } while (0)

#define ASSERT_ROUNDTRIP_1ULP(name, x, fwd_factor, inv_factor)               \
    do {                                                                      \
        double _x = (x);                                                      \
        double _y = _x * (fwd_factor);                                        \
        double _xhat = _y * (inv_factor);                                     \
        double _ulp = nextafter(fabs(_x), INFINITY) - fabs(_x);              \
        double _err = fabs(_xhat - _x);                                       \
        if (_err <= _ulp) {                                                   \
            tests_passed++;                                                   \
        } else {                                                              \
            fprintf(stderr, "FAIL %s: round-trip error %.17g > 1 ULP (%.17g)" \
                    " for x=%.17g\n", name, _err, _ulp, _x);                 \
            tests_failed++;                                                   \
        }                                                                     \
        tests_run++;                                                          \
    } while (0)

/* AC-1: Bohr constant matches CODATA 2022 */
static void test_bohr_to_nm_value(void)
{
    /* ThDD:CODATA2022:a0 — NIST CODATA 2022: a_0 = 0.052917721083(18) nm */
    ASSERT_DOUBLE_EQ("test_bohr_to_nm_value",
                     GRODFTB_BOHR_TO_NM, 0.052917721083);
}

/* AC-2: Hartree constant matches CODATA 2022 */
static void test_hartree_to_kjmol_value(void)
{
    /* ThDD:CODATA2022:Eh — E_h * N_A * 1e-3 = 2625.4996394799 kJ/mol */
    ASSERT_DOUBLE_EQ("test_hartree_to_kjmol_value",
                     GRODFTB_HARTREE_TO_KJMOL, 2625.4996394799);
}

/* AC-3: Force constant = Hartree_to_kJmol / Bohr_to_nm */
static void test_force_constant_derived(void)
{
    ASSERT_DOUBLE_EQ("test_force_constant_derived",
                     GRODFTB_FORCE_AU_TO_GMX,
                     GRODFTB_HARTREE_TO_KJMOL / GRODFTB_BOHR_TO_NM);
}

/* AC-4: Round-trip length < 1 ULP */
static void test_roundtrip_length(void)
{
    double test_values[] = {1.0, 18.897259886, 0.052917721083, 100.0, 1e-10};
    int n = sizeof(test_values) / sizeof(test_values[0]);
    for (int i = 0; i < n; i++) {
        char name[64];
        snprintf(name, sizeof(name), "test_roundtrip_length[%.6g]",
                 test_values[i]);
        ASSERT_ROUNDTRIP_1ULP(name, test_values[i],
                              GRODFTB_BOHR_TO_NM, GRODFTB_NM_TO_BOHR);
    }
}

/* AC-5: Round-trip energy < 1 ULP */
static void test_roundtrip_energy(void)
{
    double test_values[] = {1.0, -8.16, 0.001, 1000.0};
    int n = sizeof(test_values) / sizeof(test_values[0]);
    for (int i = 0; i < n; i++) {
        char name[64];
        snprintf(name, sizeof(name), "test_roundtrip_energy[%.6g]",
                 test_values[i]);
        ASSERT_ROUNDTRIP_1ULP(name, test_values[i],
                              GRODFTB_HARTREE_TO_KJMOL,
                              GRODFTB_KJMOL_TO_HARTREE);
    }
}

/* AC-6: Round-trip force < 1 ULP */
static void test_roundtrip_force(void)
{
    double test_values[] = {0.01, -0.05, 1.0, 1e-6};
    int n = sizeof(test_values) / sizeof(test_values[0]);
    for (int i = 0; i < n; i++) {
        char name[64];
        snprintf(name, sizeof(name), "test_roundtrip_force[%.6g]",
                 test_values[i]);
        ASSERT_ROUNDTRIP_1ULP(name, test_values[i],
                              GRODFTB_FORCE_AU_TO_GMX,
                              GRODFTB_FORCE_GMX_TO_AU);
    }
}

/* AC-7: coords_nm_to_bohr converts correctly */
static void test_coords_nm_to_bohr(void)
{
    double coords[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    double expected[6];
    for (int i = 0; i < 6; i++)
        expected[i] = coords[i] * GRODFTB_NM_TO_BOHR;

    grodftb_coords_nm_to_bohr(2, coords);

    int ok = 1;
    for (int i = 0; i < 6; i++) {
        if (coords[i] != expected[i]) {
            ok = 0;
            fprintf(stderr, "FAIL test_coords_nm_to_bohr: "
                    "coords[%d]=%.17g expected %.17g\n",
                    i, coords[i], expected[i]);
        }
    }
    tests_run++;
    if (ok) tests_passed++; else tests_failed++;
}

/* AC-8: coords_bohr_to_nm converts correctly */
static void test_coords_bohr_to_nm(void)
{
    double coords[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double expected[6];
    for (int i = 0; i < 6; i++)
        expected[i] = coords[i] * GRODFTB_BOHR_TO_NM;

    grodftb_coords_bohr_to_nm(2, coords);

    int ok = 1;
    for (int i = 0; i < 6; i++) {
        if (coords[i] != expected[i]) {
            ok = 0;
            fprintf(stderr, "FAIL test_coords_bohr_to_nm: "
                    "coords[%d]=%.17g expected %.17g\n",
                    i, coords[i], expected[i]);
        }
    }
    tests_run++;
    if (ok) tests_passed++; else tests_failed++;
}

/* AC-9: forces_au_to_gmx converts correctly */
static void test_forces_au_to_gmx(void)
{
    double forces[] = {0.01, -0.02, 0.03};
    double expected[3];
    for (int i = 0; i < 3; i++)
        expected[i] = forces[i] * GRODFTB_FORCE_AU_TO_GMX;

    grodftb_forces_au_to_gmx(1, forces);

    int ok = 1;
    for (int i = 0; i < 3; i++) {
        if (forces[i] != expected[i]) {
            ok = 0;
            fprintf(stderr, "FAIL test_forces_au_to_gmx: "
                    "forces[%d]=%.17g expected %.17g\n",
                    i, forces[i], expected[i]);
        }
    }
    tests_run++;
    if (ok) tests_passed++; else tests_failed++;
}

/* AC-10: forces_gmx_to_au converts correctly */
static void test_forces_gmx_to_au(void)
{
    double forces[] = {100.0, -200.0, 300.0};
    double expected[3];
    for (int i = 0; i < 3; i++)
        expected[i] = forces[i] * GRODFTB_FORCE_GMX_TO_AU;

    grodftb_forces_gmx_to_au(1, forces);

    int ok = 1;
    for (int i = 0; i < 3; i++) {
        if (forces[i] != expected[i]) {
            ok = 0;
            fprintf(stderr, "FAIL test_forces_gmx_to_au: "
                    "forces[%d]=%.17g expected %.17g\n",
                    i, forces[i], expected[i]);
        }
    }
    tests_run++;
    if (ok) tests_passed++; else tests_failed++;
}

/* AC-11: GROMACS constant note present — compile-time string check */
static void test_gromacs_constant_note(void)
{
    /*
     * Verify the header contains the GROMACS CODATA 2018 note.
     * We check that the GRODFTB_GROMACS_CODATA_NOTE macro is defined
     * (it's a string literal documenting the dual-constant convention).
     */
    const char *note = GRODFTB_GROMACS_CODATA_NOTE;
    int ok = (note != NULL && strstr(note, "GROMACS") != NULL
              && strstr(note, "CODATA 2018") != NULL);
    tests_run++;
    if (ok) {
        tests_passed++;
    } else {
        fprintf(stderr, "FAIL test_gromacs_constant_note: "
                "GRODFTB_GROMACS_CODATA_NOTE missing or wrong\n");
        tests_failed++;
    }
}

int main(void)
{
    test_bohr_to_nm_value();
    test_hartree_to_kjmol_value();
    test_force_constant_derived();
    test_roundtrip_length();
    test_roundtrip_energy();
    test_roundtrip_force();
    test_coords_nm_to_bohr();
    test_coords_bohr_to_nm();
    test_forces_au_to_gmx();
    test_forces_gmx_to_au();
    test_gromacs_constant_note();

    printf("test_units: %d run, %d passed, %d failed\n",
           tests_run, tests_passed, tests_failed);

    return tests_failed > 0 ? 1 : 0;
}
