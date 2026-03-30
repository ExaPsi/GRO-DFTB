/*
 * test_dftb3_energy.c — Validate DFTB3/3OB through DFTB+ C API
 *
 * JCTC Reviewer 2 Major Concern #7: Validate DFTB3 compatibility.
 * Compares C API energy against standalone DFTB+ for B1 water dimer
 * with DFTB3/3OB parameters (ThirdOrderFull = Yes).
 *
 * KNOWN LIMITATION: dftbp_get_gradients() causes Fortran ERROR STOP
 * with ThirdOrderFull in DFTB+ 25.1 C API. Energy and charges work;
 * gradients do not. DFTB3 MD is not supported until this is fixed
 * in a future DFTB+ release. This is a DFTB+ API limitation, not
 * a GRO-DFTB bug.
 *
 * Reference: standalone DFTB+ 25.1 (commit fd31d873)
 *   Energy = -8.1286981909 Ha (Mermin free energy)
 *   Charges: O1=-0.3410, H2=+0.1807, H3=+0.1654,
 *            O4=-0.3399, H5=+0.1676, H6=+0.1672
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <dftbplus.h>

/* Reference from standalone DFTB+ 25.1 with DFTB3/3OB on B1 water dimer */
static const double REF_ENERGY_HA = -8.1286981909;
static const double ENERGY_TOL    = 1.0e-8;  /* Hartree — bitwise match expected */

static int test_dftb3_energy(void)
{
    DftbPlus calc;
    DftbPlusInput input;
    int pass = 1;

    /* chdir to data directory (DFTB+ resolves paths from CWD) */
    char oldcwd[4096];
    getcwd(oldcwd, sizeof(oldcwd));
    if (chdir("tests/data/b1_dftb3") != 0) {
        fprintf(stderr, "FAIL: chdir to tests/data/b1_dftb3 failed\n");
        return 0;
    }

    /* Initialize DFTB+ with DFTB3/3OB HSD */
    dftbp_init(&calc, "/dev/null");
    dftbp_get_input_from_file(&calc, "dftb_in.hsd", &input);
    dftbp_process_input(&calc, &input);
    dftbp_input_final(&input);

    /* Set coordinates (same as geo.gen, converted to Bohr) */
    const double ANG2BOHR = 1.8897259886;
    double coords[18] = {
        -0.702196*ANG2BOHR, -0.056060*ANG2BOHR,  0.009942*ANG2BOHR,
        -1.022193*ANG2BOHR,  0.846776*ANG2BOHR, -0.011489*ANG2BOHR,
         0.257521*ANG2BOHR,  0.042121*ANG2BOHR,  0.005219*ANG2BOHR,
         2.220871*ANG2BOHR,  0.026717*ANG2BOHR,  0.000620*ANG2BOHR,
         2.597493*ANG2BOHR, -0.411663*ANG2BOHR,  0.766745*ANG2BOHR,
         2.593135*ANG2BOHR, -0.449496*ANG2BOHR, -0.744782*ANG2BOHR,
    };
    dftbp_set_coords(&calc, coords);

    /* Get energy (triggers SCC solve) */
    double energy;
    dftbp_get_energy(&calc, &energy);

    double energy_err = fabs(energy - REF_ENERGY_HA);
    printf("DFTB3/3OB Water Dimer (B1)\n");
    printf("  Energy (API):     %.10f Ha\n", energy);
    printf("  Energy (ref):     %.10f Ha\n", REF_ENERGY_HA);
    printf("  |error|:          %.2e Ha (tol: %.0e)\n", energy_err, ENERGY_TOL);

    if (energy_err > ENERGY_TOL) {
        printf("  Energy: FAIL\n");
        pass = 0;
    } else {
        printf("  Energy: PASS (matches standalone DFTB+)\n");
    }

    /* Get Mulliken charges */
    double charges[6];
    dftbp_get_gross_charges(&calc, charges);
    printf("  Mulliken charges:");
    double qtot = 0;
    for (int i = 0; i < 6; i++) {
        printf(" %.4f", charges[i]);
        qtot += charges[i];
    }
    printf("\n  Total charge:     %.2e (should be ~0)\n", qtot);

    /* Note on gradient limitation */
    printf("\n  Known limitation: dftbp_get_gradients() causes Fortran\n");
    printf("  ERROR STOP with ThirdOrderFull in DFTB+ 25.1 C API.\n");
    printf("  DFTB3 energy and charges work; gradients (and thus MD)\n");
    printf("  require a future DFTB+ fix. This is not a GRO-DFTB bug.\n");

    dftbp_final(&calc);
    chdir(oldcwd);
    return pass;
}

int main(void)
{
    printf("=== DFTB3/3OB Validation (Reviewer 2, MC#7) ===\n\n");

    int pass = test_dftb3_energy();

    printf("\n%s\n", pass ? "ALL TESTS PASSED" : "SOME TESTS FAILED");
    return pass ? 0 : 1;
}
