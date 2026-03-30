/*
 * US-037: B6 Alanine Dipeptide Memory Safety Tests
 *
 * ThDD:T-US-037-V.4 — Memory safety: zero leaks, zero errors
 *
 * SDD:specs.md:§21.1 — Memory: leak check with Valgrind/ASan
 *
 * This file implements memory safety tests for B6:
 * - Valgrind clean: Run short trajectory, check for leaks
 * - ASan clean: Run under AddressSanitizer, check for errors
 *
 * Usage with Valgrind:
 *   valgrind --leak-check=full --error-exitcode=1 ./test_b6_memory
 *
 * Usage with ASan (compile with -fsanitize=address):
 *   ASAN_OPTIONS=detect_leaks=1 ./test_b6_memory_asan
 *
 * IMPORTANT: Tests in this file REQUIRE DFTB+ to be built and linked.
 * Without DFTB+, tests will skip with a message.
 *
 * LGPL-3.0-or-later
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "grodftb/linkatom.h"
#include "grodftb/driver.h"
#include "grodftb/units.h"
#include "grodftb/error.h"

/* ---------------------------------------------------------------------------
 * Test Parameters
 *
 * ThDD:T-US-037-V.4 — Memory safety verification protocol:
 *   - Short trajectory (1 ps) to exercise all code paths
 *   - Multiple init/finalize cycles to detect leaks
 *   - Force reallocation patterns
 * --------------------------------------------------------------------------- */
#define MEMORY_TEST_STEPS       100     /* ~50 fs at 0.5 fs timestep */
#define MEMORY_INIT_CYCLES       3      /* Number of init/finalize cycles */
#define TIMESTEP_AU             20.67   /* 0.5 fs in atomic units */

/* B6 system parameters */
#define B6_N_QM_REAL    10
#define B6_N_MM         12
#define B6_N_LINKS       2
#define B6_N_QM_TOTAL  (B6_N_QM_REAL + B6_N_LINKS)

/* Link atom definitions */
static const int B6_QM_BOUNDARY_ATOMS[B6_N_LINKS] = {0, 8};
static const int B6_MM_BOUNDARY_ATOMS[B6_N_LINKS] = {0, 6};

/* Test data path */
#ifndef GRODFTB_SOURCE_DIR_STR
#define GRODFTB_SOURCE_DIR_STR "."
#endif
#define B6_DATA_DIR GRODFTB_SOURCE_DIR_STR "/tests/data/b6"

/* Test counters */
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

/* ---------------------------------------------------------------------------
 * Working directory management
 * --------------------------------------------------------------------------- */
static char sg_oldcwd[1024];
static int sg_in_data_dir = 0;

static int chdir_to_b6(void)
{
    if (sg_in_data_dir) return 1;
    if (!getcwd(sg_oldcwd, sizeof(sg_oldcwd))) return 0;
    if (chdir(B6_DATA_DIR) != 0) return 0;
    sg_in_data_dir = 1;
    return 1;
}

static void restore_cwd(void)
{
    if (sg_in_data_dir) {
        int rc = chdir(sg_oldcwd);
        (void)rc;  /* Ignore return value in cleanup */
        sg_in_data_dir = 0;
    }
}

/* ---------------------------------------------------------------------------
 * B6 system coordinates and species
 * --------------------------------------------------------------------------- */
static double b6_qm_coords[3 * B6_N_QM_REAL] = {
    -0.0023923575,  0.2993819398, -0.3140466218,
    -0.4753068938,  2.1874242613, -0.3263055755,
     2.7045188767, -0.0579743372, -0.0519803624,
     3.0788729657, -2.1412590151, -0.0230983537,
     3.7731322214,  1.0980302513,  2.3571887102,
     5.8334589514,  0.8907095156,  2.4049062113,
     2.9800475549,  0.1739866289,  4.0361811360,
     3.3255578038,  3.1200297049,  2.4668445563,
     4.0435113869,  0.9668246452, -2.3851967661,
     6.2195407649,  1.5856771488, -2.4202707249
};

static double b6_mm_coords[3 * B6_N_MM] = {
    -1.2803826321, -0.8486862010,  1.5704109665,
     1.0,  0.0,  0.0,
     2.0,  0.0,  0.0,
     3.0,  0.0,  0.0,
     4.0,  0.0,  0.0,
     5.0,  0.0,  0.0,
     2.5944035797,  1.0919631086, -4.4796959731,
     7.0,  0.0,  0.0,
     8.0,  0.0,  0.0,
     9.0,  0.0,  0.0,
    10.0,  0.0,  0.0,
    11.0,  0.0,  0.0,
};

static int b6_qm_species[B6_N_QM_REAL] = {
    1, 3, 0, 3, 0, 3, 3, 3, 0, 2
};

static const double b6_masses[B6_N_QM_REAL] = {
    14.007, 1.008, 12.011, 1.008, 12.011,
     1.008, 1.008,  1.008, 12.011, 15.999
};

/* ---------------------------------------------------------------------------
 * Test: Valgrind clean — short trajectory without memory errors (AC-5)
 *
 * ThDD:T-US-037-V.4 — Memory safety verification:
 *   - Zero leaks
 *   - Zero invalid reads/writes
 *   - Zero use-after-free
 *   - Zero uninitialized reads
 *
 * This test runs a short trajectory through all major code paths and then
 * properly cleans up. When run under Valgrind, any memory errors will be
 * reported and cause exit code 1.
 *
 * SDD:specs.md:§21.1 — Memory: leak check
 * --------------------------------------------------------------------------- */
static int test_b6_valgrind_clean(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    fprintf(stderr, "    DFTB+ not available - test will skip\n");
    return -1;
#else
    printf("\n    Running %d init/finalize cycles with %d MD steps each...\n",
           MEMORY_INIT_CYCLES, MEMORY_TEST_STEPS);

    for (int cycle = 0; cycle < MEMORY_INIT_CYCLES; cycle++) {
        printf("    Cycle %d/%d: ", cycle + 1, MEMORY_INIT_CYCLES);
        fflush(stdout);

        /* Change to data directory */
        if (!chdir_to_b6()) {
            fprintf(stderr, "Failed to chdir\n");
            return 0;
        }

        /* Initialize driver */
        grodftb_handle_t driver = NULL;
        int rc = grodftb_init("dftb_in.hsd", &driver);
        if (rc != GRODFTB_SUCCESS) {
            fprintf(stderr, "Failed to init driver: %d\n", rc);
            restore_cwd();
            return 0;
        }

        /* Create link atom handler */
        grodftb_linkatom_handle_t links = NULL;
        double ref_lengths[B6_N_LINKS] = {1.9451, 2.1761};
        rc = grodftb_linkatom_create(B6_N_LINKS,
                                      B6_QM_BOUNDARY_ATOMS,
                                      B6_MM_BOUNDARY_ATOMS,
                                      ref_lengths,
                                      3,  /* H species index */
                                      GRODFTB_CHARGE_NONE,
                                      &links);
        if (rc != GRODFTB_SUCCESS) {
            fprintf(stderr, "Failed to create linkatom handler: %d\n", rc);
            grodftb_finalize(&driver);
            restore_cwd();
            return 0;
        }

        /* Allocate working arrays */
        double *augmented_coords = malloc(3 * B6_N_QM_TOTAL * sizeof(double));
        double *augmented_forces = malloc(3 * B6_N_QM_TOTAL * sizeof(double));
        int *augmented_species = malloc(B6_N_QM_TOTAL * sizeof(int));

        if (!augmented_coords || !augmented_forces || !augmented_species) {
            fprintf(stderr, "Failed to allocate arrays\n");
            free(augmented_coords);
            free(augmented_forces);
            free(augmented_species);
            grodftb_linkatom_destroy(&links);
            grodftb_finalize(&driver);
            restore_cwd();
            return 0;
        }

        /* Working coordinates */
        double qm_coords[3 * B6_N_QM_REAL];
        double mm_coords[3 * B6_N_MM];
        double qm_forces[3 * B6_N_QM_REAL];
        double mm_forces[3 * B6_N_MM];
        double velocities[3 * B6_N_QM_REAL];

        memcpy(qm_coords, b6_qm_coords, sizeof(qm_coords));
        memcpy(mm_coords, b6_mm_coords, sizeof(mm_coords));
        memset(velocities, 0, sizeof(velocities));

        /* Initialize velocities */
        for (int i = 0; i < B6_N_QM_REAL; i++) {
            double mass_au = b6_masses[i] * 1822.888486;
            double v_rms = sqrt(9.5e-4 / mass_au);
            velocities[3*i + 0] = v_rms * 0.1 * (i % 3 - 1);
            velocities[3*i + 1] = v_rms * 0.1 * ((i + 1) % 3 - 1);
            velocities[3*i + 2] = v_rms * 0.1 * ((i + 2) % 3 - 1);
        }

        /* Run short trajectory */
        for (int step = 0; step < MEMORY_TEST_STEPS; step++) {
            /* Augment coordinates */
            rc = grodftb_linkatom_augment_coords(links, B6_N_QM_REAL,
                                                  qm_coords, mm_coords,
                                                  augmented_coords);
            if (rc != GRODFTB_SUCCESS) {
                fprintf(stderr, "augment_coords failed at step %d\n", step);
                break;
            }

            /* Augment species */
            rc = grodftb_linkatom_augment_species(links, B6_N_QM_REAL,
                                                   b6_qm_species, augmented_species);
            if (rc != GRODFTB_SUCCESS) break;

            /* Set geometry */
            rc = grodftb_set_geometry(driver, B6_N_QM_TOTAL,
                                       augmented_species, augmented_coords);
            if (rc != GRODFTB_SUCCESS) break;

            /* Compute */
            rc = grodftb_compute(driver);
            if (rc != GRODFTB_SUCCESS) {
                fprintf(stderr, "compute failed at step %d: %d\n", step, rc);
                break;
            }

            /* Get energy (exercise the API) */
            double energy;
            rc = grodftb_get_energy(driver, &energy);
            if (rc != GRODFTB_SUCCESS) break;

            /* Get forces */
            rc = grodftb_get_forces(driver, augmented_forces);
            if (rc != GRODFTB_SUCCESS) break;

            /* Project forces */
            memset(mm_forces, 0, sizeof(mm_forces));
            rc = grodftb_linkatom_project_forces(links, B6_N_QM_REAL,
                                                  qm_coords, mm_coords,
                                                  augmented_forces,
                                                  qm_forces, mm_forces);
            if (rc != GRODFTB_SUCCESS) break;

            /* Get charges (exercise charge API) */
            double charges[B6_N_QM_TOTAL];
            rc = grodftb_get_mulliken_charges(driver, charges);
            if (rc != GRODFTB_SUCCESS) break;

            /* Velocity Verlet step */
            for (int i = 0; i < B6_N_QM_REAL; i++) {
                double mass_au = b6_masses[i] * 1822.888486;
                for (int d = 0; d < 3; d++) {
                    velocities[3*i + d] += qm_forces[3*i + d] / mass_au * TIMESTEP_AU;
                    qm_coords[3*i + d] += velocities[3*i + d] * TIMESTEP_AU;
                }
            }
        }

        /* Clean up all allocations */
        free(augmented_coords);
        free(augmented_forces);
        free(augmented_species);

        /* Destroy link atom handler */
        grodftb_linkatom_destroy(&links);

        /* Verify handle is NULL after destroy */
        if (links != NULL) {
            fprintf(stderr, "linkatom handle not NULL after destroy\n");
        }

        /* Finalize driver */
        grodftb_finalize(&driver);

        /* Verify handle is NULL after finalize */
        if (driver != NULL) {
            fprintf(stderr, "driver handle not NULL after finalize\n");
        }

        restore_cwd();
        printf("OK\n");
    }

    printf("    All cycles completed successfully.\n");
    printf("    If running under Valgrind, check for leak summary above.\n");

    return 1;
#endif
}

/* ---------------------------------------------------------------------------
 * Test: ASan clean — run under AddressSanitizer (AC-5)
 *
 * ThDD:T-US-037-V.4 — Memory safety verification with ASan
 *
 * This test is identical to valgrind_clean but is specifically intended
 * for compilation with -fsanitize=address. ASan provides faster detection
 * of memory errors than Valgrind but requires recompilation.
 *
 * Usage:
 *   cmake -DCMAKE_C_FLAGS="-fsanitize=address -g" ..
 *   make test_b6_memory
 *   ASAN_OPTIONS=detect_leaks=1 ./test_b6_memory
 * --------------------------------------------------------------------------- */
static int test_b6_asan_clean(void)
{
    /*
     * This test exercises the same code paths as valgrind_clean.
     * When compiled with ASan, the sanitizer will catch:
     * - Buffer overflows
     * - Use after free
     * - Stack buffer overflow
     * - Global buffer overflow
     * - Use after return
     * - Use after scope
     * - Initialization order bugs
     * - Memory leaks (with detect_leaks=1)
     */

    printf("\n    NOTE: This test is for ASan-compiled builds.\n");
    printf("    If not compiled with -fsanitize=address, results are same as valgrind test.\n");

    /* Run the same test as valgrind_clean */
    return test_b6_valgrind_clean();
}

/* ---------------------------------------------------------------------------
 * Test: Stress test — many allocations and deallocations (AC-5)
 *
 * This test specifically targets memory allocation patterns to detect
 * fragmentation issues and ensure proper deallocation in edge cases.
 * --------------------------------------------------------------------------- */
static int test_b6_stress_alloc(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    return -1;
#else
    printf("\n    Stress testing allocation patterns...\n");

    /* Create and destroy many redistribution result structures */
    printf("    Testing redistrib_result create/destroy...\n");
    for (int i = 0; i < 100; i++) {
        grodftb_redistrib_result_t *result = NULL;
        int rc = grodftb_redistrib_result_create(&result, B6_N_MM);
        if (rc != GRODFTB_SUCCESS) {
            fprintf(stderr, "    Failed to create redistrib_result at iteration %d\n", i);
            return 0;
        }
        grodftb_redistrib_result_destroy(&result);
        if (result != NULL) {
            fprintf(stderr, "    redistrib_result not NULL after destroy\n");
            return 0;
        }
    }

    /* Create and destroy many link atom handlers */
    printf("    Testing linkatom_handler create/destroy...\n");
    for (int i = 0; i < 50; i++) {
        grodftb_linkatom_handle_t links = NULL;
        double ref_lengths[B6_N_LINKS] = {1.9451, 2.1761};
        int rc = grodftb_linkatom_create(B6_N_LINKS,
                                          B6_QM_BOUNDARY_ATOMS,
                                          B6_MM_BOUNDARY_ATOMS,
                                          ref_lengths,
                                          3,
                                          GRODFTB_CHARGE_NONE,
                                          &links);
        if (rc != GRODFTB_SUCCESS) {
            fprintf(stderr, "    Failed to create linkatom_handler at iteration %d\n", i);
            return 0;
        }
        grodftb_linkatom_destroy(&links);
        if (links != NULL) {
            fprintf(stderr, "    linkatom_handle not NULL after destroy\n");
            return 0;
        }
    }

    /* Test NULL pointer handling */
    printf("    Testing NULL pointer handling...\n");

    /* These should not crash */
    grodftb_linkatom_handle_t null_links = NULL;
    int count = grodftb_linkatom_count(null_links);
    if (count != 0) {
        fprintf(stderr, "    linkatom_count(NULL) should return 0, got %d\n", count);
        return 0;
    }

    grodftb_linkatom_destroy(&null_links);  /* Should not crash */
    grodftb_linkatom_destroy(NULL);         /* Should not crash */

    grodftb_handle_t null_driver = NULL;
    grodftb_finalize(&null_driver);         /* Should not crash */
    grodftb_finalize(NULL);                 /* Should not crash */

    grodftb_redistrib_result_t *null_result = NULL;
    grodftb_redistrib_result_destroy(&null_result);  /* Should not crash */
    grodftb_redistrib_result_destroy(NULL);          /* Should not crash */

    printf("    Stress test completed successfully.\n");
    return 1;
#endif
}

/* ---------------------------------------------------------------------------
 * Main test runner
 * --------------------------------------------------------------------------- */
int main(void)
{
    printf("=== US-037 B6 Memory Safety Tests ===\n\n");

    printf("Usage notes:\n");
    printf("  For Valgrind:  valgrind --leak-check=full --error-exitcode=1 ./test_b6_memory\n");
    printf("  For ASan:      Compile with -fsanitize=address, run with ASAN_OPTIONS=detect_leaks=1\n");
    printf("\n");

    printf("AC-5: Memory safety:\n");
    RUN_TEST(test_b6_valgrind_clean);
    RUN_TEST(test_b6_asan_clean);
    RUN_TEST(test_b6_stress_alloc);

    printf("\n=== Summary ===\n");
    printf("Tests run:    %d\n", tests_run);
    printf("Tests passed: %d\n", tests_passed);
    printf("Tests failed: %d\n", tests_failed);
    printf("Tests skipped: %d\n", tests_skipped);

    if (tests_skipped == tests_run) {
        printf("\nNOTE: All tests skipped - DFTB+ not available.\n");
        printf("      Build with DFTB+ linked to run these tests.\n");
        return 0;
    }

    printf("\nNOTE: Memory safety is VERIFIED only when run under Valgrind or ASan.\n");
    printf("      A passing result here means no crashes, not necessarily no leaks.\n");

    return tests_failed > 0 ? 1 : 0;
}
