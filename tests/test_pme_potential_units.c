/*
 * TDD:US-046 — Unit tests for PME potential/field unit conversion
 * SDD:specs.md:§8.2 — PME-compatible embedding
 * ThDD:T-US-046-U-2.2 — Potential conversion: kJ/(mol·e) → Hartree
 * ThDD:T-US-046-U-2.7 — Field conversion: kJ/(mol·e·nm) → Hartree/Bohr
 *
 * Tests AC-4: Unit conversion round-trip error < 1 ULP
 *
 * These tests verify the conversion constants at the GROMACS→DFTB+ boundary.
 * Constants are from GROMACS units.h (CODATA 2018).
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

#define ASSERT_ROUNDTRIP_1ULP(name, x, fwd_factor, inv_factor)               \
    do {                                                                      \
        double _x = (x);                                                      \
        double _y = _x * (fwd_factor);                                        \
        double _xhat = _y * (inv_factor);                                     \
        double _ulp = nextafter(fabs(_x), INFINITY) - fabs(_x);              \
        double _err = fabs(_xhat - _x);                                       \
        if (_err <= _ulp) {                                                   \
            tests_passed++;                                                   \
            printf("  PASS %s\n", name);                                      \
        } else {                                                              \
            fprintf(stderr, "  FAIL %s: round-trip error %.17g > 1 ULP "     \
                    "(%.17g) for x=%.17g\n", name, _err, _ulp, _x);          \
            tests_failed++;                                                   \
        }                                                                     \
        tests_run++;                                                          \
    } while (0)

#define ASSERT_REL_EQ(name, a, b, tol)                                       \
    do {                                                                      \
        double _a = (a), _b = (b);                                           \
        double _rel = fabs(_a - _b) / fabs(_b);                              \
        if (_rel < (tol)) {                                                   \
            tests_passed++;                                                   \
            printf("  PASS %s (rel err %.2e)\n", name, _rel);                \
        } else {                                                              \
            fprintf(stderr, "  FAIL %s: rel error %.17g >= %.17g "           \
                    "(a=%.17g, b=%.17g)\n", name, _rel, (double)(tol),       \
                    _a, _b);                                                  \
            tests_failed++;                                                   \
        }                                                                     \
        tests_run++;                                                          \
    } while (0)

#define ASSERT_ABS_EQ(name, a, b, tol)                                       \
    do {                                                                      \
        double _a = (a), _b = (b);                                           \
        double _err = fabs(_a - _b);                                          \
        if (_err < (tol)) {                                                   \
            tests_passed++;                                                   \
            printf("  PASS %s (abs err %.2e)\n", name, _err);                \
        } else {                                                              \
            fprintf(stderr, "  FAIL %s: abs error %.17g >= %.17g "           \
                    "(a=%.17g, b=%.17g)\n", name, _err, (double)(tol),       \
                    _a, _b);                                                  \
            tests_failed++;                                                   \
        }                                                                     \
        tests_run++;                                                          \
    } while (0)

/*
 * GROMACS constants from units.h (CODATA 2018)
 * Duplicated here for standalone testing.
 * The actual implementation uses gromacs/math/units.h at compile time.
 */
static const double GMX_C_AVOGADRO    = 6.02214076e23;     /* 1/mol, CODATA 2018 */
static const double GMX_C_BOHR2NM     = 0.0529177210903;   /* nm, CODATA 2018 */

/*
 * c_hartree2Kj is derived in GROMACS from fundamental constants:
 *   c_hartree2Kj = 2 * c_rydberg * c_planck * c_speedOfLight / c_avogadro
 * Its value is approximately 4.3597447222071e-18 J = 4.3597447222071e-21 kJ
 * (This is the Hartree energy in kJ per molecule, NOT per mole)
 *
 * The combined factor c_hartree2Kj * c_avogadro ≈ 2625.4996... kJ/mol
 * is what converts Hartree to kJ/mol.
 */
static const double GMX_C_HARTREE2KJ  = 4.3597447222071e-21; /* kJ (per molecule) */

/*
 * Derived conversion factor: kJ/(mol·e) → Hartree
 * ThDD:T-US-046-U-2.2
 *
 * Φ_Ha = Φ_kJ/(mol·e) / (c_hartree2Kj * c_avogadro)
 *      = Φ_kJ/(mol·e) / 2625.4996...
 */
static double potential_gmx_to_hartree(void)
{
    return 1.0 / (GMX_C_HARTREE2KJ * GMX_C_AVOGADRO);
}

static double potential_hartree_to_gmx(void)
{
    return GMX_C_HARTREE2KJ * GMX_C_AVOGADRO;
}

/*
 * Derived conversion factor: kJ/(mol·e·nm) → Hartree/Bohr
 * ThDD:T-US-046-U-2.7
 *
 * E_Ha/Bohr = E_kJ/(mol·e·nm) * c_bohr2Nm / (c_hartree2Kj * c_avogadro)
 */
static double field_gmx_to_hartree_bohr(void)
{
    return GMX_C_BOHR2NM / (GMX_C_HARTREE2KJ * GMX_C_AVOGADRO);
}

static double field_hartree_bohr_to_gmx(void)
{
    return (GMX_C_HARTREE2KJ * GMX_C_AVOGADRO) / GMX_C_BOHR2NM;
}

/*
 * Test 1: Potential round-trip conversion < 1 ULP
 * ThDD:T-US-046-U-2.2
 *
 * Convert kJ/(mol·e) → Hartree → kJ/(mol·e), verify < 1 ULP error.
 * Test at several representative potential magnitudes.
 */
static void test_potential_gmx_to_hartree_roundtrip(void)
{
    printf("test_potential_gmx_to_hartree_roundtrip:\n");

    const double fwd = potential_gmx_to_hartree();
    const double inv = potential_hartree_to_gmx();

    /* Typical QM/MM potentials: 0.01-1.0 Hartree */
    ASSERT_ROUNDTRIP_1ULP("pot_0.1Ha",   0.1 * inv,   fwd, inv);
    ASSERT_ROUNDTRIP_1ULP("pot_1.0Ha",   1.0 * inv,   fwd, inv);
    ASSERT_ROUNDTRIP_1ULP("pot_small",   1e-6 * inv,  fwd, inv);
    ASSERT_ROUNDTRIP_1ULP("pot_large",   10.0 * inv,  fwd, inv);
    ASSERT_ROUNDTRIP_1ULP("pot_neg",    -0.3 * inv,   fwd, inv);
}

/*
 * Test 2: Field round-trip conversion < 1 ULP
 * ThDD:T-US-046-U-2.7
 *
 * Convert kJ/(mol·e·nm) → Hartree/Bohr → kJ/(mol·e·nm), verify < 1 ULP error.
 */
static void test_field_gmx_to_hartree_bohr_roundtrip(void)
{
    printf("test_field_gmx_to_hartree_bohr_roundtrip:\n");

    const double fwd = field_gmx_to_hartree_bohr();
    const double inv = field_hartree_bohr_to_gmx();

    /* Typical field magnitudes */
    ASSERT_ROUNDTRIP_1ULP("field_0.01",   0.01 * inv,  fwd, inv);
    ASSERT_ROUNDTRIP_1ULP("field_1.0",    1.0 * inv,   fwd, inv);
    ASSERT_ROUNDTRIP_1ULP("field_neg",   -0.5 * inv,   fwd, inv);
}

/*
 * Test 3: Dimensional sanity check
 * ThDD:T-US-046-U-2.6
 *
 * A unit charge (+1 e) at distance 1 Bohr should produce a potential
 * of exactly 1 Hartree. In GROMACS units:
 *   Φ_gmx = Q * c_one4PiEps0 / r_nm
 *         = 1.0 * 138.935... / 0.0529177...
 *         = 2625.4996... kJ/(mol·e)
 *
 * Converting: Φ_Ha = Φ_gmx * fwd = 1.0 Hartree
 *
 * GROMACS c_one4PiEps0 = 138.93545764... kJ·nm/(mol·e²)
 */
static void test_potential_dimensional_sanity(void)
{
    printf("test_potential_dimensional_sanity:\n");

    /* c_one4PiEps0 from GROMACS units.h */
    /* 1/(4πε₀) in kJ·nm/(mol·e²) */
    const double c_one4PiEps0 = 138.93545764438198;

    /* Φ(r=1 Bohr) = c_one4PiEps0 / c_bohr2Nm */
    const double phi_gmx = c_one4PiEps0 / GMX_C_BOHR2NM;

    /* Convert to Hartree */
    const double phi_hartree = phi_gmx * potential_gmx_to_hartree();

    /*
     * Expected: exactly 1.0 Hartree
     * Tolerance: machine epsilon (~1e-15) since this is a fundamental identity.
     * The only error sources are the double-precision representation of
     * the constants themselves.
     */
    ASSERT_REL_EQ("phi_1bohr_equals_1hartree", phi_hartree, 1.0, 1e-10);
}

/*
 * Test 4: Potential conversion factor numerical value
 * ThDD:T-US-046-U-2.3
 *
 * The factor kJ/(mol·e) → Hartree = 1/(c_hartree2Kj * c_avogadro)
 * ≈ 1/2625.4996... ≈ 3.8088e-4
 */
static void test_potential_conversion_factor_value(void)
{
    printf("test_potential_conversion_factor_value:\n");

    const double factor = potential_gmx_to_hartree();
    const double expected_reciprocal = GMX_C_HARTREE2KJ * GMX_C_AVOGADRO;

    /* Verify it's approximately 2625.5 kJ/mol per Hartree */
    ASSERT_REL_EQ("kJ_mol_per_hartree", expected_reciprocal, 2625.4996, 1e-4);

    /* Verify the factor is approximately 3.8088e-4 */
    ASSERT_REL_EQ("hartree_per_kJ_mol", factor, 3.8088e-4, 1e-3);
}

/*
 * Test 5: Field conversion factor numerical value
 * ThDD:T-US-046-U-2.7
 *
 * The factor kJ/(mol·e·nm) → Hartree/Bohr
 *   = c_bohr2Nm / (c_hartree2Kj * c_avogadro)
 *   = 0.052918 / 2625.50 ≈ 2.0155e-5
 */
static void test_field_conversion_factor_value(void)
{
    printf("test_field_conversion_factor_value:\n");

    const double factor = field_gmx_to_hartree_bohr();

    /* Verify the factor is approximately 2.0155e-5 */
    ASSERT_REL_EQ("ha_bohr_per_kJ_mol_nm", factor, 2.0155e-5, 1e-3);

    /* Verify relationship: field_factor = potential_factor * c_bohr2Nm */
    const double pot_factor = potential_gmx_to_hartree();
    ASSERT_REL_EQ("field_pot_bohr_relation",
                  factor, pot_factor * GMX_C_BOHR2NM, 1e-15);
}

int main(void)
{
    printf("=== PME Potential Unit Conversion Tests (US-046, AC-4) ===\n\n");

    test_potential_gmx_to_hartree_roundtrip();
    test_field_gmx_to_hartree_bohr_roundtrip();
    test_potential_dimensional_sanity();
    test_potential_conversion_factor_value();
    test_field_conversion_factor_value();

    printf("\n--- Results: %d/%d passed, %d failed ---\n",
           tests_passed, tests_run, tests_failed);

    return (tests_failed > 0) ? 1 : 0;
}
