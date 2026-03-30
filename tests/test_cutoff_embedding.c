/*
 * ============================================================================
 * DEPRECATED: This file is retained for historical documentation only.
 * ============================================================================
 *
 * This test file was written during the TDD Red Phase for US-038 to test
 * a hypothetical standalone libgrodftb API for cutoff embedding. However,
 * the actual implementation was done in the GROMACS backend instead:
 *
 *   - DFTBForceProvider::gatherEmbeddingCharges()
 *   - DFTBForceProvider::computeEmbeddingMmForces()
 *
 * The embedding functionality is validated by:
 *   1. tests/test_embedding_fd.c -- FD validation through libgrodftb API
 *   2. external/gromacs/src/gromacs/applied_forces/qmmm/tests/
 *      dftbforceprovider_embedding.cpp -- GROMACS integration tests
 *
 * The tests in this file call stub functions (grodftb_gather_embedding_atoms,
 * grodftb_compute_embedding_potential) that were never implemented and
 * always return GRODFTB_ERR_NOT_INITIALIZED.
 *
 * DO NOT add this file to the build. It is kept for documentation of the
 * original test design intent during the TDD Red Phase.
 *
 * Date deprecated: 2026-02-04
 * Reason: API design changed during US-038 implementation
 * See: docs/reports/US-038/report.md for full verification report
 * ============================================================================
 *
 * Original header:
 *
 * SDD:specs.md:S8 -- TDD Red Phase tests for US-038: Cutoff-Based Electrostatic Embedding
 *
 * US-038: Cutoff embedding in DFTBForceProvider::calculateForces()
 *
 * ThDD:T-US-038-N.1 -- Cutoff neighbor gathering algorithm
 * ThDD:T-US-038-1.1 -- MM electrostatic potential at QM sites
 * ThDD:T-US-038-V.1 -- Energy agreement with exact Coulomb sum
 * ThDD:T-US-038-V.5 -- No QM-QM double-counting
 *
 * These tests are written BEFORE implementation (TDD Red Phase).
 * All tests should initially FAIL because the cutoff embedding module
 * does not exist yet.
 *
 * Reference data provenance:
 *   - B4 benchmark: tests/data/b4/ (water dimer + 1 point charge)
 *   - DFTB+ version: 25.1 (commit fd31d873)
 *   - SK parameters: mio-1-1
 *   - Total energy: -8.18640458887598 Ha (from B4 reference.json)
 *   - MM forces: (0.011319515597, -0.006946588945, 0.000215086973) Ha/Bohr
 *
 * LGPL-3.0-or-later
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "grodftb/driver.h"
#include "grodftb/error.h"
#include "grodftb/units.h"

/* ---------------------------------------------------------------------------
 * Test framework macros and counters
 * --------------------------------------------------------------------------- */
static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

#define RUN_TEST(test_func) do { \
    printf("  Running %s... ", #test_func); \
    fflush(stdout); \
    tests_run++; \
    if (test_func()) { \
        printf("PASS\n"); \
        tests_passed++; \
    } else { \
        printf("FAIL\n"); \
        tests_failed++; \
    } \
} while(0)

/* ---------------------------------------------------------------------------
 * ThDD:T-US-038-N.1 -- Cutoff neighbor gathering parameters
 *
 * Default cutoff: 1.2 nm = 22.677 Bohr
 * Unit conversions use GRODFTB_NM_TO_BOHR from units.h
 * --------------------------------------------------------------------------- */
#define DEFAULT_CUTOFF_NM       1.2
#define DEFAULT_CUTOFF_BOHR     (DEFAULT_CUTOFF_NM * GRODFTB_NM_TO_BOHR)

/* ---------------------------------------------------------------------------
 * ThDD:T-US-038-V.1 -- Energy tolerance from specs.md S21.1
 * SDD:specs.md:S21.1 -- Energy match < 10^-8 Hartree
 * --------------------------------------------------------------------------- */
#define ENERGY_TOLERANCE_HA     1e-8

/* ---------------------------------------------------------------------------
 * Reference data from B4 benchmark (tests/data/b4/reference.json)
 *
 * Provenance: DFTB+ 25.1, commit fd31d873, mio-1-1 SK parameters
 * Generated: 2026-02-01T03:01:19Z by tools/generate_reference_data.sh
 *
 * ThDD:T-US-038-V.1 -- Expected energy for embedded water dimer
 * --------------------------------------------------------------------------- */
#define B4_ENERGY_HA           (-8.18640458887598)

/* B4 MM charge (single point charge at -3.0 Angstrom from origin along x) */
#define B4_MM_CHARGE            1.0
#define B4_MM_POS_X_ANGSTROM  (-3.0)
#define B4_MM_POS_Y_ANGSTROM    0.0
#define B4_MM_POS_Z_ANGSTROM    0.0

/* Convert Angstrom to Bohr: 1 Angstrom = 1.8897259886 Bohr */
#define ANGSTROM_TO_BOHR        1.8897259886

/* B4 MM forces on external charge (from detailed.out line 89)
 * ThDD:T-US-038-2.2 -- Coulomb back-reaction force formula */
#define B4_MM_FORCE_X_HA_BOHR   0.011319515597
#define B4_MM_FORCE_Y_HA_BOHR (-0.006946588945)
#define B4_MM_FORCE_Z_HA_BOHR   0.000215086973

/* ---------------------------------------------------------------------------
 * TDD Red Phase: Stub implementations for unimplemented API
 *
 * These stubs allow the tests to compile and run during the Red Phase.
 * They return GRODFTB_ERR_NOT_INITIALIZED to indicate "not implemented yet".
 *
 * When the real implementation is added to libgrodftb, these stubs will
 * be replaced by the library symbols during linking.
 *
 * SDD:specs.md:S8.1 -- Cutoff embedding algorithm interface
 * --------------------------------------------------------------------------- */

/*
 * Gather MM atoms within cutoff radius of QM region.
 *
 * ThDD:T-US-038-N.1 -- Cutoff neighbor gathering:
 *   For each MM atom J, include if min_A(|R_A - R_J|) < r_cutoff
 *   where A ranges over all QM atoms.
 *
 * @param n_qm            Number of QM atoms
 * @param qm_coords       QM atom coordinates in Bohr, flat [3*n_qm]
 * @param n_mm            Number of MM atoms
 * @param mm_coords       MM atom coordinates in Bohr, flat [3*n_mm]
 * @param cutoff_bohr     Cutoff radius in Bohr
 * @param n_gathered_out  Output: number of MM atoms within cutoff
 * @param mm_indices_out  Output: indices of gathered MM atoms [n_gathered_out]
 *                        Caller must pre-allocate with size >= n_mm
 * @return GRODFTB_SUCCESS on success
 *
 * NOTE: This function signature is for TDD Red Phase. The actual
 * implementation may differ (e.g., in the GROMACS backend).
 */
#ifndef GRODFTB_HAS_GATHER_EMBEDDING
static int grodftb_gather_embedding_atoms(
    int n_qm, const double *qm_coords,
    int n_mm, const double *mm_coords,
    double cutoff_bohr,
    int *n_gathered_out, int *mm_indices_out)
{
    /* TDD Red Phase stub: function not implemented yet */
    (void)n_qm; (void)qm_coords; (void)n_mm; (void)mm_coords;
    (void)cutoff_bohr; (void)n_gathered_out; (void)mm_indices_out;
    return GRODFTB_ERR_NOT_INITIALIZED;
}
#endif

/*
 * Compute electrostatic potential at QM sites from MM charges.
 *
 * ThDD:T-US-038-1.1 -- MM potential formula:
 *   Phi_MM(R_A) = sum_{J in MM} Q_J / |R_A - R_J|
 *
 * @param n_qm        Number of QM atoms
 * @param qm_coords   QM coordinates in Bohr, flat [3*n_qm]
 * @param n_mm        Number of MM atoms
 * @param mm_charges  MM charges in elementary charge [n_mm]
 * @param mm_coords   MM coordinates in Bohr, flat [3*n_mm]
 * @param potential   Output: potential at each QM site in Hartree [n_qm]
 * @return GRODFTB_SUCCESS on success
 */
#ifndef GRODFTB_HAS_COMPUTE_POTENTIAL
static int grodftb_compute_embedding_potential(
    int n_qm, const double *qm_coords,
    int n_mm, const double *mm_charges, const double *mm_coords,
    double *potential)
{
    /* TDD Red Phase stub: function not implemented yet */
    (void)n_qm; (void)qm_coords; (void)n_mm; (void)mm_charges;
    (void)mm_coords; (void)potential;
    return GRODFTB_ERR_NOT_INITIALIZED;
}
#endif

/* ---------------------------------------------------------------------------
 * Test helper: Path resolution for test data
 * --------------------------------------------------------------------------- */
static const char *get_b4_hsd_path(void)
{
    static char path[1024];
    const char *srcdir = getenv("GRODFTB_SOURCE_DIR");
    if (!srcdir) {
        srcdir = GRODFTB_SOURCE_DIR_STR;
    }
    snprintf(path, sizeof(path), "%s/tests/data/b4/dftb_in.hsd", srcdir);
    return path;
}

static const char *get_b4_dir(void)
{
    static char path[1024];
    const char *srcdir = getenv("GRODFTB_SOURCE_DIR");
    if (!srcdir) {
        srcdir = GRODFTB_SOURCE_DIR_STR;
    }
    snprintf(path, sizeof(path), "%s/tests/data/b4", srcdir);
    return path;
}

/* ---------------------------------------------------------------------------
 * Test helper: Distance calculation
 * --------------------------------------------------------------------------- */
static double distance_3d(const double *a, const double *b)
{
    double dx = a[0] - b[0];
    double dy = a[1] - b[1];
    double dz = a[2] - b[2];
    return sqrt(dx*dx + dy*dy + dz*dz);
}

static double min_qm_mm_distance(int n_qm, const double *qm_coords,
                                  const double *mm_pos)
{
    double min_dist = 1e30;
    for (int i = 0; i < n_qm; i++) {
        double d = distance_3d(&qm_coords[3*i], mm_pos);
        if (d < min_dist) min_dist = d;
    }
    return min_dist;
}

/* ---------------------------------------------------------------------------
 * Test 1: test_cutoff_gathers_within_radius
 *
 * ThDD:T-US-038-N.1 -- Cutoff neighbor gathering algorithm
 * SDD:US-038.md:AC-1 -- MM charges within cutoff gathered correctly
 *
 * Set up: 1 QM atom at origin, 3 MM atoms at distances 0.5, 1.0, 1.5 nm
 * Cutoff = 1.2 nm
 * Assert: MM atoms at 0.5 and 1.0 nm are in embedding list, 1.5 nm is not
 * --------------------------------------------------------------------------- */
static int test_cutoff_gathers_within_radius(void)
{
    /* QM atom at origin */
    double qm_coords[3] = {0.0, 0.0, 0.0};
    int n_qm = 1;

    /* MM atoms at 0.5, 1.0, 1.5 nm along x axis */
    double mm_coords[9] = {
        0.5 * GRODFTB_NM_TO_BOHR, 0.0, 0.0,   /* 0.5 nm - should be gathered */
        1.0 * GRODFTB_NM_TO_BOHR, 0.0, 0.0,   /* 1.0 nm - should be gathered */
        1.5 * GRODFTB_NM_TO_BOHR, 0.0, 0.0    /* 1.5 nm - should NOT be gathered */
    };
    int n_mm = 3;

    double cutoff_bohr = 1.2 * GRODFTB_NM_TO_BOHR;

    int n_gathered = 0;
    int mm_indices[3] = {-1, -1, -1};

    /*
     * RED PHASE: This call should fail because grodftb_gather_embedding_atoms
     * is not implemented yet.
     */
    int rc = grodftb_gather_embedding_atoms(
        n_qm, qm_coords,
        n_mm, mm_coords,
        cutoff_bohr,
        &n_gathered, mm_indices
    );

    if (rc != GRODFTB_SUCCESS) {
        printf("\n    grodftb_gather_embedding_atoms returned error %d\n", rc);
        return 0;
    }

    /* Verify exactly 2 atoms gathered (indices 0 and 1) */
    if (n_gathered != 2) {
        printf("\n    Expected 2 atoms gathered, got %d\n", n_gathered);
        return 0;
    }

    /* Verify the correct indices (order may vary) */
    int found_0 = 0, found_1 = 0, found_2 = 0;
    for (int i = 0; i < n_gathered; i++) {
        if (mm_indices[i] == 0) found_0 = 1;
        if (mm_indices[i] == 1) found_1 = 1;
        if (mm_indices[i] == 2) found_2 = 1;
    }

    if (!found_0 || !found_1) {
        printf("\n    Expected atoms 0 and 1 in gathered list\n");
        return 0;
    }
    if (found_2) {
        printf("\n    Atom 2 (1.5 nm) should NOT be in gathered list\n");
        return 0;
    }

    return 1;
}

/* ---------------------------------------------------------------------------
 * Test 2: test_cutoff_excludes_beyond_radius
 *
 * ThDD:T-US-038-N.1 -- Cutoff neighbor gathering algorithm
 * SDD:US-038.md:AC-2 -- MM charges beyond cutoff excluded
 *
 * Set up: 1 QM atom, MM atoms at 1.3, 1.5, 2.0 nm
 * Cutoff = 1.2 nm
 * Assert: No MM atoms in embedding list
 * --------------------------------------------------------------------------- */
static int test_cutoff_excludes_beyond_radius(void)
{
    /* QM atom at origin */
    double qm_coords[3] = {0.0, 0.0, 0.0};
    int n_qm = 1;

    /* MM atoms all beyond cutoff */
    double mm_coords[9] = {
        1.3 * GRODFTB_NM_TO_BOHR, 0.0, 0.0,
        1.5 * GRODFTB_NM_TO_BOHR, 0.0, 0.0,
        2.0 * GRODFTB_NM_TO_BOHR, 0.0, 0.0
    };
    int n_mm = 3;

    double cutoff_bohr = 1.2 * GRODFTB_NM_TO_BOHR;

    int n_gathered = -1;
    int mm_indices[3];

    int rc = grodftb_gather_embedding_atoms(
        n_qm, qm_coords,
        n_mm, mm_coords,
        cutoff_bohr,
        &n_gathered, mm_indices
    );

    if (rc != GRODFTB_SUCCESS) {
        printf("\n    grodftb_gather_embedding_atoms returned error %d\n", rc);
        return 0;
    }

    if (n_gathered != 0) {
        printf("\n    Expected 0 atoms gathered, got %d\n", n_gathered);
        return 0;
    }

    return 1;
}

/* ---------------------------------------------------------------------------
 * Test 3: test_cutoff_at_boundary
 *
 * ThDD:T-US-038-N.1 -- Cutoff uses strict < comparison
 *
 * MM atom exactly at cutoff distance should be EXCLUDED (d < cutoff, not <=)
 * --------------------------------------------------------------------------- */
static int test_cutoff_at_boundary(void)
{
    /* QM atom at origin */
    double qm_coords[3] = {0.0, 0.0, 0.0};
    int n_qm = 1;

    /* MM atom exactly at cutoff distance */
    double cutoff_bohr = 1.2 * GRODFTB_NM_TO_BOHR;
    double mm_coords[3] = {cutoff_bohr, 0.0, 0.0};
    int n_mm = 1;

    int n_gathered = -1;
    int mm_indices[1];

    int rc = grodftb_gather_embedding_atoms(
        n_qm, qm_coords,
        n_mm, mm_coords,
        cutoff_bohr,
        &n_gathered, mm_indices
    );

    if (rc != GRODFTB_SUCCESS) {
        printf("\n    grodftb_gather_embedding_atoms returned error %d\n", rc);
        return 0;
    }

    /* ThDD:T-US-038-N.1 -- strict < comparison means exactly at cutoff is excluded */
    if (n_gathered != 0) {
        printf("\n    Atom at exactly cutoff should be excluded (strict < comparison)\n");
        printf("    Expected 0 gathered, got %d\n", n_gathered);
        return 0;
    }

    return 1;
}

/* ---------------------------------------------------------------------------
 * Test 4: test_potential_b4_single_charge
 *
 * ThDD:T-US-038-1.1 -- Potential at QM sites matches libgrodftb direct computation
 * SDD:US-038.md:AC-3 -- Potential matches to < 10^-12 relative error
 *
 * B4 benchmark: water dimer + 1 point charge
 * --------------------------------------------------------------------------- */
static int test_potential_b4_single_charge(void)
{
    /*
     * B4 water dimer QM geometry (from tests/data/b4/geo.gen, in Angstrom)
     * Convert to Bohr for internal calculations
     *
     * ThDD:06_theory.md:Eq1.1 -- Coulomb potential formula:
     *   Phi(R_A) = sum_J Q_J / |R_A - R_J|
     */

    /* QM atom coordinates from B4 geo.gen (converted to Bohr) */
    double qm_coords_ang[18] = {
        -0.702196, -0.056060,  0.009942,  /* O1 */
        -1.022193,  0.846776, -0.011489,  /* H1 */
         0.257521,  0.042121,  0.005219,  /* H2 */
         2.220871,  0.026717,  0.000620,  /* O2 */
         2.597493, -0.411663,  0.766745,  /* H3 */
         2.593135, -0.449496, -0.744782   /* H4 */
    };
    double qm_coords[18];
    for (int i = 0; i < 18; i++) {
        qm_coords[i] = qm_coords_ang[i] * ANGSTROM_TO_BOHR;
    }
    int n_qm = 6;

    /* MM point charge from B4 dftb_in.hsd: Q=1.0 at (-3.0, 0.0, 0.0) Angstrom */
    double mm_charges[1] = {B4_MM_CHARGE};
    double mm_coords[3] = {
        B4_MM_POS_X_ANGSTROM * ANGSTROM_TO_BOHR,
        B4_MM_POS_Y_ANGSTROM * ANGSTROM_TO_BOHR,
        B4_MM_POS_Z_ANGSTROM * ANGSTROM_TO_BOHR
    };
    int n_mm = 1;

    /* Output potential array */
    double potential[6];

    /*
     * RED PHASE: This call should fail because grodftb_compute_embedding_potential
     * is not implemented yet.
     */
    int rc = grodftb_compute_embedding_potential(
        n_qm, qm_coords,
        n_mm, mm_charges, mm_coords,
        potential
    );

    if (rc != GRODFTB_SUCCESS) {
        printf("\n    grodftb_compute_embedding_potential returned error %d\n", rc);
        return 0;
    }

    /*
     * Verify potential by direct calculation:
     * Phi_A = Q_J / |R_A - R_J|
     *
     * ThDD:T-US-038-1.1 -- Validate computed potential against analytic formula
     */
    int all_match = 1;
    double max_rel_error = 0.0;

    for (int i = 0; i < n_qm; i++) {
        double r = distance_3d(&qm_coords[3*i], mm_coords);
        double expected = mm_charges[0] / r;

        double rel_error = fabs((potential[i] - expected) / expected);
        if (rel_error > max_rel_error) max_rel_error = rel_error;

        if (rel_error > 1e-12) {
            printf("\n    Atom %d: potential = %.15e, expected = %.15e, rel_error = %.2e\n",
                   i, potential[i], expected, rel_error);
            all_match = 0;
        }
    }

    if (!all_match) {
        printf("    Max relative error: %.2e (tolerance: 1e-12)\n", max_rel_error);
        return 0;
    }

    return 1;
}

/* ---------------------------------------------------------------------------
 * Test 5: test_energy_b4_embedded
 *
 * ThDD:T-US-038-V.1 -- Energy matches reference to < 10^-8 Ha
 * SDD:specs.md:S21.1 -- Energy match < 10^-8 Hartree
 *
 * B4 benchmark: water dimer + 1 point charge
 * Expected: -8.18640458887598 Ha (from tests/data/b4/reference.json)
 * --------------------------------------------------------------------------- */
static int test_energy_b4_embedded(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    printf(" SKIPPED (no DFTB+)");
    return 1;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b4_hsd_path();
    char oldcwd[1024];
    const char *datadir = get_b4_dir();

    /* DFTB+ requires running from the data directory (relative paths in HSD) */
    if (!getcwd(oldcwd, sizeof(oldcwd))) {
        printf("\n    getcwd failed\n");
        return 0;
    }
    if (chdir(datadir) != 0) {
        printf("\n    chdir to b4 data dir failed\n");
        return 0;
    }

    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS) {
        chdir(oldcwd);
        printf("\n    grodftb_init failed: %d (%s)\n", rc, grodftb_error_string(rc));
        return 0;
    }

    /* Compute (geometry and embedding are set in the HSD for B4) */
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        chdir(oldcwd);
        printf("\n    grodftb_compute failed: %d\n", rc);
        return 0;
    }

    /* Get energy */
    double energy = 0.0;
    rc = grodftb_get_energy(handle, &energy);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        chdir(oldcwd);
        printf("\n    grodftb_get_energy failed: %d\n", rc);
        return 0;
    }

    grodftb_finalize(&handle);
    chdir(oldcwd);

    /*
     * ThDD:T-US-038-V.1 -- Compare against archived reference
     *
     * Reference: B4_ENERGY_HA = -8.18640458887598 Ha
     * Provenance: tests/data/b4/reference.json
     * DFTB+ version: 25.1 (commit fd31d873)
     * SK parameters: mio-1-1
     */
    double error = fabs(energy - B4_ENERGY_HA);
    if (error > ENERGY_TOLERANCE_HA) {
        printf("\n    Energy mismatch:\n");
        printf("      computed: %.15e Ha\n", energy);
        printf("      expected: %.15e Ha\n", B4_ENERGY_HA);
        printf("      error:    %.6e Ha (tolerance: %.0e)\n", error, ENERGY_TOLERANCE_HA);
        return 0;
    }

    printf(" (error: %.2e Ha)", error);
    return 1;
#endif
}

/* ---------------------------------------------------------------------------
 * Test 6: test_no_qm_in_embedding
 *
 * ThDD:T-US-038-V.5 -- Verify QM indices never appear in embedding list
 *
 * The embedding list should contain only MM atoms. QM-QM Coulomb interactions
 * are handled internally by DFTB+ (gamma_AB matrix).
 * --------------------------------------------------------------------------- */
static int test_no_qm_in_embedding(void)
{
    /*
     * Set up a system where QM and MM coordinates overlap at some positions.
     * The gathering function must exclude positions that match QM atoms.
     *
     * This test simulates a scenario where the caller might accidentally
     * include QM atom positions in the MM list.
     */

    /* QM atoms */
    double qm_coords[6] = {
        0.0, 0.0, 0.0,    /* QM atom 0 */
        2.0, 0.0, 0.0     /* QM atom 1 */
    };
    int n_qm = 2;

    /*
     * MM atoms: one at a different position, one at exactly the QM position.
     * The one at the QM position should be excluded (it's a QM atom).
     */
    double mm_coords[6] = {
        1.0, 0.0, 0.0,    /* MM atom 0: valid, 1 Bohr from QM 0 */
        0.0, 0.0, 0.0     /* MM atom 1: invalid, same as QM atom 0 */
    };
    int n_mm = 2;

    /* QM indices that should be excluded from embedding */
    int qm_indices[2] = {0, 1};

    double cutoff_bohr = 5.0;  /* Large enough to include everything geometrically */

    int n_gathered = -1;
    int mm_indices[2];

    /*
     * The gathering function should have an additional parameter to exclude
     * specific indices. For now, test that it at least handles the case.
     *
     * NOTE: The actual API may differ. This test documents the expected behavior.
     */
    int rc = grodftb_gather_embedding_atoms(
        n_qm, qm_coords,
        n_mm, mm_coords,
        cutoff_bohr,
        &n_gathered, mm_indices
    );

    if (rc != GRODFTB_SUCCESS) {
        printf("\n    grodftb_gather_embedding_atoms returned error %d\n", rc);
        return 0;
    }

    /*
     * ThDD:T-US-038-V.5 -- No double-counting
     *
     * If the implementation detects overlapping positions and excludes them,
     * only MM atom 0 should be gathered.
     *
     * Alternative: The API may require the caller to ensure no overlap.
     * In that case, this test documents the expected caller contract.
     */
    printf(" (API-dependent: checking for overlap detection)");

    /* For now, just verify the call succeeded. Full verification requires
     * the actual implementation to define the overlap-handling behavior. */
    return 1;
}

/* ---------------------------------------------------------------------------
 * Test 7: test_empty_embedding_region
 *
 * ThDD:T-US-038-V.5 -- Gas-phase behavior when no MM atoms within cutoff
 *
 * When no MM atoms are within the cutoff, the embedding energy should be
 * zero and the system should behave as gas-phase.
 * --------------------------------------------------------------------------- */
static int test_empty_embedding_region(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    printf(" SKIPPED (no DFTB+)");
    return 1;
#else
    /*
     * Test that when we call grodftb_set_embedding_charges with ncharges=0,
     * we get gas-phase results.
     *
     * Reference: B1 water dimer gas-phase energy
     * (We compare embedded with ncharges=0 against reference gas-phase)
     */

    const char *srcdir = getenv("GRODFTB_SOURCE_DIR");
    if (!srcdir) srcdir = GRODFTB_SOURCE_DIR_STR;

    char hsd_path[1024];
    snprintf(hsd_path, sizeof(hsd_path), "%s/tests/data/b1/dftb_in.hsd", srcdir);

    char datadir[1024];
    snprintf(datadir, sizeof(datadir), "%s/tests/data/b1", srcdir);

    char oldcwd[1024];
    if (!getcwd(oldcwd, sizeof(oldcwd))) {
        printf("\n    getcwd failed\n");
        return 0;
    }
    if (chdir(datadir) != 0) {
        printf("\n    chdir to b1 data dir failed\n");
        return 0;
    }

    grodftb_handle_t handle = NULL;
    int rc = grodftb_init(hsd_path, &handle);
    if (rc != GRODFTB_SUCCESS) {
        chdir(oldcwd);
        printf("\n    grodftb_init failed: %d\n", rc);
        return 0;
    }

    /* Set embedding with no charges (gas-phase) */
    rc = grodftb_set_embedding_charges(handle, 0, NULL, NULL);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        chdir(oldcwd);
        printf("\n    grodftb_set_embedding_charges(0, NULL, NULL) failed: %d\n", rc);
        return 0;
    }

    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        chdir(oldcwd);
        printf("\n    grodftb_compute failed: %d\n", rc);
        return 0;
    }

    double energy = 0.0;
    rc = grodftb_get_energy(handle, &energy);
    grodftb_finalize(&handle);
    chdir(oldcwd);

    if (rc != GRODFTB_SUCCESS) {
        printf("\n    grodftb_get_energy failed: %d\n", rc);
        return 0;
    }

    /*
     * B1 reference energy (gas-phase water dimer)
     * From tests/data/b1/reference.json: -8.16019724583451 Ha
     *
     * Provenance: DFTB+ 25.1, commit fd31d873, mio-1-1 SK parameters
     */
    double b1_ref_energy = -8.16019724583451;
    double error = fabs(energy - b1_ref_energy);

    if (error > ENERGY_TOLERANCE_HA) {
        printf("\n    Gas-phase energy mismatch:\n");
        printf("      computed: %.15e Ha\n", energy);
        printf("      expected: %.15e Ha\n", b1_ref_energy);
        printf("      error:    %.6e Ha\n", error);
        return 0;
    }

    printf(" (gas-phase energy matches B1 reference)");
    return 1;
#endif
}

/* ---------------------------------------------------------------------------
 * Main test runner
 * --------------------------------------------------------------------------- */
int main(void)
{
    printf("=== US-038 Cutoff Embedding Tests (TDD Red Phase) ===\n\n");

    printf("Cutoff neighbor gathering (ThDD:T-US-038-N.1):\n");
    RUN_TEST(test_cutoff_gathers_within_radius);
    RUN_TEST(test_cutoff_excludes_beyond_radius);
    RUN_TEST(test_cutoff_at_boundary);

    printf("\nPotential computation (ThDD:T-US-038-1.1):\n");
    RUN_TEST(test_potential_b4_single_charge);

    printf("\nEnergy validation (ThDD:T-US-038-V.1):\n");
    RUN_TEST(test_energy_b4_embedded);

    printf("\nSanity checks (ThDD:T-US-038-V.5):\n");
    RUN_TEST(test_no_qm_in_embedding);
    RUN_TEST(test_empty_embedding_region);

    printf("\n=== Summary ===\n");
    printf("Tests run:    %d\n", tests_run);
    printf("Tests passed: %d\n", tests_passed);
    printf("Tests failed: %d\n", tests_failed);

    if (tests_failed > 0) {
        printf("\nNOTE: This is the TDD Red Phase. Tests are expected to fail\n");
        printf("      until the cutoff embedding module is implemented.\n");
    }

    return tests_failed > 0 ? 1 : 0;
}
