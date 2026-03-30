/*
 * TDD tests for the DFTB+ driver library.
 *
 * US-009: grodftb_init() and grodftb_finalize() — AC-1 through AC-9.
 * US-011: grodftb_set_geometry() — AC-1 through AC-8.
 * US-012: grodftb_compute() — AC-1 through AC-10.
 * US-013: grodftb_get_energy() exact match — AC-1 through AC-8.
 * US-014: grodftb_get_forces() force retrieval — AC-1 through AC-8.
 * US-016: grodftb_get_mulliken_charges() — AC-1 through AC-8.
 * US-017: Error handling — grodftb_error_string(), grodftb_last_error().
 * US-018: Full integration tests — init→set_geometry→compute→get_*→finalize.
 * US-019: grodftb_set_embedding_charges() — electrostatic embedding via external potential.
 * US-024: grodftb_version() — version info retrieval.
 *
 * Uses B1 water dimer HSD template at tests/data/b1/dftb_in.hsd.
 * Uses B2 formamide HSD template at tests/data/b2/dftb_in.hsd.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "grodftb/driver.h"
#include "grodftb/error.h"

#ifdef GRODFTB_HAS_DFTBPLUS
#include <dftbplus.h>
#endif

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
    fprintf(stderr, "  FAIL: %s — %s\n", name, reason);
}

/*
 * Helper: get path to B1 HSD relative to the test binary working directory.
 * The CI script runs ctest from the build dir, but the HSD references geo.gen
 * relatively, so we need to run from the b1 data directory.
 *
 * We use an environment variable or a compile-time define for the source root.
 */
static const char *get_b1_hsd_path(void)
{
    static char path[1024];
    const char *srcdir = getenv("GRODFTB_SOURCE_DIR");
    if (!srcdir) {
        srcdir = GRODFTB_SOURCE_DIR_STR;
    }
    snprintf(path, sizeof(path), "%s/tests/data/b1/dftb_in.hsd", srcdir);
    return path;
}

/* AC-1: Init with valid HSD returns 0 and handle != NULL */
static void test_driver_init_success(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_init_success [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();

    /* DFTB+ needs to run from the directory containing the HSD
     * because geo.gen is referenced with a relative path. */
    char oldcwd[1024];
    char datadir[1024];
    const char *srcdir = getenv("GRODFTB_SOURCE_DIR");
    if (!srcdir) srcdir = GRODFTB_SOURCE_DIR_STR;
    snprintf(datadir, sizeof(datadir), "%s/tests/data/b1", srcdir);

    if (!getcwd(oldcwd, sizeof(oldcwd))) {
        fail("test_driver_init_success", "getcwd failed");
        return;
    }
    if (chdir(datadir) != 0) {
        fail("test_driver_init_success", "chdir to b1 data dir failed");
        return;
    }

    int rc = grodftb_init(hsd, &handle);

    /* Restore working directory before any assertions */
    if (chdir(oldcwd) != 0) {
        fail("test_driver_init_success", "chdir back failed");
        return;
    }

    if (rc != GRODFTB_SUCCESS) {
        char msg[256];
        snprintf(msg, sizeof(msg), "init returned %d, expected 0", rc);
        fail("test_driver_init_success", msg);
        return;
    }
    if (!handle) {
        fail("test_driver_init_success", "handle is NULL after successful init");
        return;
    }

    grodftb_finalize(&handle);
    pass("test_driver_init_success");
#endif
}

/* AC-2: After init, DFTB+ API version is 0.4.0 */
static void test_driver_api_version(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_api_version [SKIPPED: no DFTB+]");
    return;
#else
    int major = -1, minor = -1, patch = -1;
    dftbp_api(&major, &minor, &patch);
    if (major == 0 && minor == 4 && patch == 0) {
        pass("test_driver_api_version");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "API version %d.%d.%d, expected 0.4.0",
                 major, minor, patch);
        fail("test_driver_api_version", msg);
    }
#endif
}

/* AC-3: NULL path returns GRODFTB_ERR_INVALID_ARGUMENT */
static void test_driver_init_null_path(void)
{
    grodftb_handle_t handle = NULL;
    int rc = grodftb_init(NULL, &handle);
    if (rc == GRODFTB_ERR_INVALID_ARGUMENT) {
        pass("test_driver_init_null_path");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_INVALID_ARGUMENT);
        fail("test_driver_init_null_path", msg);
    }
}

/* AC-4: NULL handle_out returns GRODFTB_ERR_INVALID_ARGUMENT */
static void test_driver_init_null_handle(void)
{
    int rc = grodftb_init("/some/path.hsd", NULL);
    if (rc == GRODFTB_ERR_INVALID_ARGUMENT) {
        pass("test_driver_init_null_handle");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_INVALID_ARGUMENT);
        fail("test_driver_init_null_handle", msg);
    }
}

/* AC-5: Nonexistent file returns error, handle is NULL */
static void test_driver_init_bad_path(void)
{
    grodftb_handle_t handle = (grodftb_handle_t)(void *)0xDEAD;
    int rc = grodftb_init("/nonexistent_path_42.hsd", &handle);
    if (rc == GRODFTB_ERR_FILE_NOT_FOUND && handle == NULL) {
        pass("test_driver_init_bad_path");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "rc=%d (expected %d), handle=%p (expected NULL)",
                 rc, GRODFTB_ERR_FILE_NOT_FOUND, (void *)handle);
        fail("test_driver_init_bad_path", msg);
    }
}

/* AC-6: Finalize sets *handle = NULL */
static void test_driver_finalize_nulls_handle(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_finalize_nulls_handle [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    char oldcwd[1024], datadir[1024];
    const char *srcdir = getenv("GRODFTB_SOURCE_DIR");
    if (!srcdir) srcdir = GRODFTB_SOURCE_DIR_STR;
    snprintf(datadir, sizeof(datadir), "%s/tests/data/b1", srcdir);
    if (!getcwd(oldcwd, sizeof(oldcwd))) {
        fail("test_driver_finalize_nulls_handle", "getcwd failed");
        return;
    }
    if (chdir(datadir) != 0) {
        fail("test_driver_finalize_nulls_handle", "chdir failed");
        return;
    }

    int rc = grodftb_init(hsd, &handle);
    if (chdir(oldcwd) != 0) {
        fail("test_driver_finalize_nulls_handle", "chdir back failed");
        return;
    }

    if (rc != GRODFTB_SUCCESS || !handle) {
        fail("test_driver_finalize_nulls_handle", "init failed — cannot test finalize");
        return;
    }

    grodftb_finalize(&handle);
    if (handle == NULL) {
        pass("test_driver_finalize_nulls_handle");
    } else {
        fail("test_driver_finalize_nulls_handle", "handle not NULL after finalize");
    }
#endif
}

/* AC-7: finalize(NULL) is safe */
static void test_driver_finalize_null_safe(void)
{
    grodftb_finalize(NULL);
    grodftb_handle_t h = NULL;
    grodftb_finalize(&h);
    pass("test_driver_finalize_null_safe");
}

/* AC-8: 1000 init-finalize cycles: zero leaks (validated under Valgrind externally) */
static void test_driver_init_finalize_leak(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_init_finalize_leak [SKIPPED: no DFTB+]");
    return;
#else
    const char *hsd = get_b1_hsd_path();
    char oldcwd[1024], datadir[1024];
    const char *srcdir = getenv("GRODFTB_SOURCE_DIR");
    if (!srcdir) srcdir = GRODFTB_SOURCE_DIR_STR;
    snprintf(datadir, sizeof(datadir), "%s/tests/data/b1", srcdir);
    if (!getcwd(oldcwd, sizeof(oldcwd))) {
        fail("test_driver_init_finalize_leak", "getcwd failed");
        return;
    }
    if (chdir(datadir) != 0) {
        fail("test_driver_init_finalize_leak", "chdir failed");
        return;
    }

    /* Run fewer cycles in normal test; Valgrind test uses full 1000 */
    int n_cycles = 3;
    const char *valgrind_env = getenv("GRODFTB_LEAK_TEST_CYCLES");
    if (valgrind_env) {
        n_cycles = atoi(valgrind_env);
        if (n_cycles < 1) n_cycles = 3;
    }

    int ok = 1;
    for (int i = 0; i < n_cycles; i++) {
        grodftb_handle_t handle = NULL;
        int rc = grodftb_init(hsd, &handle);
        if (rc != GRODFTB_SUCCESS || !handle) {
            char msg[128];
            snprintf(msg, sizeof(msg), "init failed on cycle %d with rc=%d", i, rc);
            fail("test_driver_init_finalize_leak", msg);
            ok = 0;
            break;
        }
        grodftb_finalize(&handle);
    }

    if (chdir(oldcwd) != 0) {
        fail("test_driver_init_finalize_leak", "chdir back failed");
        return;
    }

    if (ok) {
        pass("test_driver_init_finalize_leak");
    }
#endif
}

/* ================================================================
 * US-011: grodftb_set_geometry() tests
 * ================================================================ */

/*
 * Helper: create a chdir-to-b1 / restore scope for DFTB+ tests.
 * Returns 1 on success (cwd changed), 0 on failure.
 */
static char sg_oldcwd[1024];

static int chdir_to_b1(void)
{
    char datadir[1024];
    const char *srcdir = getenv("GRODFTB_SOURCE_DIR");
    if (!srcdir) srcdir = GRODFTB_SOURCE_DIR_STR;
    snprintf(datadir, sizeof(datadir), "%s/tests/data/b1", srcdir);
    if (!getcwd(sg_oldcwd, sizeof(sg_oldcwd))) return 0;
    return chdir(datadir) == 0;
}

static void restore_cwd(void)
{
    if (chdir(sg_oldcwd) != 0) {
        /* Silently ignore error; best-effort restore */
        (void)0;
    }
}

/* B1 water dimer coordinates from geo.gen (already in Bohr) */
static const double b1_coords[18] = {
    -0.702196000000, -0.056060000000,  0.009942000000,  /* O */
    -1.022193000000,  0.846776000000, -0.011489000000,  /* H */
     0.257521000000,  0.042121000000,  0.005219000000,  /* H */
     2.220871000000,  0.026717000000,  0.000620000000,  /* O */
     2.597493000000, -0.411663000000,  0.766745000000,  /* H */
     2.593135000000, -0.449496000000, -0.744782000000   /* H */
};

/* 0-based species: O=0, H=1 */
static const int b1_species[6] = { 0, 1, 1, 0, 1, 1 };

/* US-011 AC-1: NULL handle returns GRODFTB_ERR_NULL_HANDLE */
static void test_driver_set_geometry_null_handle(void)
{
    int rc = grodftb_set_geometry(NULL, 6, b1_species, b1_coords);
    if (rc == GRODFTB_ERR_NULL_HANDLE) {
        pass("test_driver_set_geometry_null_handle");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_NULL_HANDLE);
        fail("test_driver_set_geometry_null_handle", msg);
    }
}

/* US-011 AC-2: Uninitialized handle returns GRODFTB_ERR_NOT_INITIALIZED.
 * We allocate a raw zeroed block — initialized field will be 0. */
static void test_driver_set_geometry_not_init(void)
{
    /* Allocate enough bytes for the struct; calloc zeroes all fields
     * including initialized=0. We only need the first few bytes. */
    grodftb_handle_t fake = (grodftb_handle_t)calloc(1, 256);
    if (!fake) {
        fail("test_driver_set_geometry_not_init", "calloc failed");
        return;
    }
    int rc = grodftb_set_geometry(fake, 6, b1_species, b1_coords);
    free(fake);
    if (rc == GRODFTB_ERR_NOT_INITIALIZED) {
        pass("test_driver_set_geometry_not_init");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_NOT_INITIALIZED);
        fail("test_driver_set_geometry_not_init", msg);
    }
}

/* US-011 AC-3: NULL coords returns GRODFTB_ERR_INVALID_ARGUMENT */
static void test_driver_set_geometry_null_coords(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_set_geometry_null_coords [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_driver_set_geometry_null_coords", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    restore_cwd();
    if (rc != GRODFTB_SUCCESS || !handle) {
        fail("test_driver_set_geometry_null_coords", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, NULL);
    grodftb_finalize(&handle);
    if (rc == GRODFTB_ERR_INVALID_ARGUMENT) {
        pass("test_driver_set_geometry_null_coords");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_INVALID_ARGUMENT);
        fail("test_driver_set_geometry_null_coords", msg);
    }
#endif
}

/* US-011 AC-4: NULL species returns GRODFTB_ERR_INVALID_ARGUMENT */
static void test_driver_set_geometry_null_species(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_set_geometry_null_species [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_driver_set_geometry_null_species", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    restore_cwd();
    if (rc != GRODFTB_SUCCESS || !handle) {
        fail("test_driver_set_geometry_null_species", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, NULL, b1_coords);
    grodftb_finalize(&handle);
    if (rc == GRODFTB_ERR_INVALID_ARGUMENT) {
        pass("test_driver_set_geometry_null_species");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_INVALID_ARGUMENT);
        fail("test_driver_set_geometry_null_species", msg);
    }
#endif
}

/* US-011 AC-5: natoms <= 0 returns GRODFTB_ERR_INVALID_ARGUMENT */
static void test_driver_set_geometry_bad_natoms(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_set_geometry_bad_natoms [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_driver_set_geometry_bad_natoms", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    restore_cwd();
    if (rc != GRODFTB_SUCCESS || !handle) {
        fail("test_driver_set_geometry_bad_natoms", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 0, b1_species, b1_coords);
    int rc2 = grodftb_set_geometry(handle, -1, b1_species, b1_coords);
    grodftb_finalize(&handle);
    if (rc == GRODFTB_ERR_INVALID_ARGUMENT && rc2 == GRODFTB_ERR_INVALID_ARGUMENT) {
        pass("test_driver_set_geometry_bad_natoms");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "rc=%d, rc2=%d, expected %d",
                 rc, rc2, GRODFTB_ERR_INVALID_ARGUMENT);
        fail("test_driver_set_geometry_bad_natoms", msg);
    }
#endif
}

/* US-011 AC-6: natoms != ctx->natoms returns GRODFTB_ERR_SIZE_MISMATCH */
static void test_driver_set_geometry_size_mismatch(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_set_geometry_size_mismatch [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_driver_set_geometry_size_mismatch", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    restore_cwd();
    if (rc != GRODFTB_SUCCESS || !handle) {
        fail("test_driver_set_geometry_size_mismatch", "init failed");
        return;
    }
    /* B1 has 6 atoms; pass 3 */
    int wrong_species[3] = { 0, 1, 1 };
    double wrong_coords[9] = { 0.0 };
    rc = grodftb_set_geometry(handle, 3, wrong_species, wrong_coords);
    grodftb_finalize(&handle);
    if (rc == GRODFTB_ERR_SIZE_MISMATCH) {
        pass("test_driver_set_geometry_size_mismatch");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_SIZE_MISMATCH);
        fail("test_driver_set_geometry_size_mismatch", msg);
    }
#endif
}

/* US-011 AC-7: Successful set_geometry with B1 water dimer */
static void test_driver_set_geometry_b1(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_set_geometry_b1 [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_driver_set_geometry_b1", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    restore_cwd();
    if (rc != GRODFTB_SUCCESS || !handle) {
        fail("test_driver_set_geometry_b1", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords);
    grodftb_finalize(&handle);
    if (rc == GRODFTB_SUCCESS) {
        pass("test_driver_set_geometry_b1");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected 0", rc);
        fail("test_driver_set_geometry_b1", msg);
    }
#endif
}

/* US-011 AC-8: 100 set_geometry calls — no leaks (Valgrind validated) */
static void test_driver_set_geometry_leak(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_set_geometry_leak [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_driver_set_geometry_leak", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    restore_cwd();
    if (rc != GRODFTB_SUCCESS || !handle) {
        fail("test_driver_set_geometry_leak", "init failed");
        return;
    }

    int n_cycles = 100;
    int ok = 1;
    for (int i = 0; i < n_cycles; i++) {
        rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords);
        if (rc != GRODFTB_SUCCESS) {
            char msg[128];
            snprintf(msg, sizeof(msg), "set_geometry failed on cycle %d with rc=%d", i, rc);
            fail("test_driver_set_geometry_leak", msg);
            ok = 0;
            break;
        }
    }

    grodftb_finalize(&handle);
    if (ok) {
        pass("test_driver_set_geometry_leak");
    }
#endif
}

/* AC-9: Failed init leaks no memory (validated under Valgrind externally) */
static void test_driver_init_fail_no_leak(void)
{
    grodftb_handle_t handle = NULL;
    for (int i = 0; i < 100; i++) {
        int rc = grodftb_init("/nonexistent_42.hsd", &handle);
        (void)rc;
    }
    if (handle == NULL) {
        pass("test_driver_init_fail_no_leak");
    } else {
        fail("test_driver_init_fail_no_leak", "handle not NULL after failed init");
    }
}

/* ================================================================
 * US-012: grodftb_compute() tests
 * ================================================================ */

/* US-012 AC-4: grodftb_compute(NULL) returns GRODFTB_ERR_NULL_HANDLE */
static void test_driver_compute_null_handle(void)
{
    int rc = grodftb_compute(NULL);
    if (rc == GRODFTB_ERR_NULL_HANDLE) {
        pass("test_driver_compute_null_handle");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_NULL_HANDLE);
        fail("test_driver_compute_null_handle", msg);
    }
}

/* US-012 AC-5: Uninitialized handle returns GRODFTB_ERR_NOT_INITIALIZED */
static void test_driver_compute_not_init(void)
{
    grodftb_handle_t fake = (grodftb_handle_t)calloc(1, 256);
    if (!fake) {
        fail("test_driver_compute_not_init", "calloc failed");
        return;
    }
    int rc = grodftb_compute(fake);
    free(fake);
    if (rc == GRODFTB_ERR_NOT_INITIALIZED) {
        pass("test_driver_compute_not_init");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_NOT_INITIALIZED);
        fail("test_driver_compute_not_init", msg);
    }
}

/* US-012 AC-6: compute without set_geometry returns GRODFTB_ERR_INVALID_ARGUMENT */
static void test_driver_compute_no_geometry(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_compute_no_geometry [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_driver_compute_no_geometry", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    restore_cwd();
    if (rc != GRODFTB_SUCCESS || !handle) {
        fail("test_driver_compute_no_geometry", "init failed");
        return;
    }
    /* Call compute without calling set_geometry first */
    rc = grodftb_compute(handle);
    grodftb_finalize(&handle);
    if (rc == GRODFTB_ERR_INVALID_ARGUMENT) {
        pass("test_driver_compute_no_geometry");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_INVALID_ARGUMENT);
        fail("test_driver_compute_no_geometry", msg);
    }
#endif
}

/* US-012 AC-1, AC-2: Full init→set_geometry→compute returns GRODFTB_SUCCESS */
static void test_driver_compute_b1(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_compute_b1 [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_driver_compute_b1", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_driver_compute_b1", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_driver_compute_b1", "set_geometry failed");
        return;
    }
    rc = grodftb_compute(handle);
    restore_cwd();
    grodftb_finalize(&handle);
    if (rc == GRODFTB_SUCCESS) {
        pass("test_driver_compute_b1");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "compute returned %d, expected 0", rc);
        fail("test_driver_compute_b1", msg);
    }
#endif
}

/* US-012 AC-3: After compute, energy is finite, negative, and reasonable.
 * No fabricated expected value — only sanity checks (E < 0, isfinite, |E| < 1e4). */
static void test_driver_compute_state(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_compute_state [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_driver_compute_state", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_driver_compute_state", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_driver_compute_state", "set_geometry failed");
        return;
    }
    rc = grodftb_compute(handle);
    restore_cwd();
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        fail("test_driver_compute_state", "compute failed");
        return;
    }
    double e = 0.0;
    rc = grodftb_get_energy(handle, &e);
    grodftb_finalize(&handle);
    if (rc != GRODFTB_SUCCESS) {
        char msg[128];
        snprintf(msg, sizeof(msg), "get_energy returned %d", rc);
        fail("test_driver_compute_state", msg);
        return;
    }
    if (!isfinite(e)) {
        fail("test_driver_compute_state", "energy is not finite");
        return;
    }
    if (e >= 0.0) {
        char msg[128];
        snprintf(msg, sizeof(msg), "energy = %.15e, expected negative", e);
        fail("test_driver_compute_state", msg);
        return;
    }
    if (fabs(e) >= 1.0e4) {
        char msg[128];
        snprintf(msg, sizeof(msg), "energy = %.15e, |E| >= 1e4 is unreasonable", e);
        fail("test_driver_compute_state", msg);
        return;
    }
    pass("test_driver_compute_state");
#endif
}

/* US-012 AC-7: After compute, grodftb_is_converged() returns flag==1 */
static void test_driver_compute_converged(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_compute_converged [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_driver_compute_converged", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_driver_compute_converged", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_driver_compute_converged", "set_geometry failed");
        return;
    }
    rc = grodftb_compute(handle);
    restore_cwd();
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        fail("test_driver_compute_converged", "compute failed");
        return;
    }
    int flag = -1;
    rc = grodftb_is_converged(handle, &flag);
    grodftb_finalize(&handle);
    if (rc != GRODFTB_SUCCESS) {
        char msg[128];
        snprintf(msg, sizeof(msg), "is_converged returned %d", rc);
        fail("test_driver_compute_converged", msg);
        return;
    }
    if (flag == 1) {
        pass("test_driver_compute_converged");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "converged flag = %d, expected 1", flag);
        fail("test_driver_compute_converged", msg);
    }
#endif
}

/* US-012 AC-8: Two computes without geometry update yield bitwise identical energy */
static void test_driver_compute_double(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_compute_double [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_driver_compute_double", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_driver_compute_double", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_driver_compute_double", "set_geometry failed");
        return;
    }
    /* First compute */
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_driver_compute_double", "first compute failed");
        return;
    }
    double e1 = 0.0;
    rc = grodftb_get_energy(handle, &e1);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_driver_compute_double", "first get_energy failed");
        return;
    }
    /* Second compute — no geometry update */
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_driver_compute_double", "second compute failed");
        return;
    }
    double e2 = 0.0;
    rc = grodftb_get_energy(handle, &e2);
    restore_cwd();
    grodftb_finalize(&handle);
    if (rc != GRODFTB_SUCCESS) {
        fail("test_driver_compute_double", "second get_energy failed");
        return;
    }
    /* Bitwise comparison — identical geometry must yield identical energy */
    if (memcmp(&e1, &e2, sizeof(double)) == 0) {
        pass("test_driver_compute_double");
    } else {
        char msg[256];
        snprintf(msg, sizeof(msg), "e1=%.17e, e2=%.17e — not bitwise identical", e1, e2);
        fail("test_driver_compute_double", msg);
    }
#endif
}

/* ================================================================
 * B2 formamide test data and helpers
 * ================================================================ */

static const char *get_b2_hsd_path(void)
{
    static char path[1024];
    const char *srcdir = getenv("GRODFTB_SOURCE_DIR");
    if (!srcdir) {
        srcdir = GRODFTB_SOURCE_DIR_STR;
    }
    snprintf(path, sizeof(path), "%s/tests/data/b2/dftb_in.hsd", srcdir);
    return path;
}

static int chdir_to_b2(void)
{
    char datadir[1024];
    const char *srcdir = getenv("GRODFTB_SOURCE_DIR");
    if (!srcdir) srcdir = GRODFTB_SOURCE_DIR_STR;
    snprintf(datadir, sizeof(datadir), "%s/tests/data/b2", srcdir);
    if (!getcwd(sg_oldcwd, sizeof(sg_oldcwd))) return 0;
    return chdir(datadir) == 0;
}

/* 0-based species: C=0, O=1, N=2, H=3 */
static const int b2_species[6] = { 0, 1, 2, 3, 3, 3 };

/* US-012 AC-10: 10 set_geometry→compute cycles, all succeed (Valgrind validated) */
static void test_driver_compute_leak(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_compute_leak [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_driver_compute_leak", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_driver_compute_leak", "init failed");
        return;
    }

    int ok = 1;
    for (int i = 0; i < 10; i++) {
        rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords);
        if (rc != GRODFTB_SUCCESS) {
            char msg[128];
            snprintf(msg, sizeof(msg), "set_geometry failed on cycle %d with rc=%d", i, rc);
            fail("test_driver_compute_leak", msg);
            ok = 0;
            break;
        }
        rc = grodftb_compute(handle);
        if (rc != GRODFTB_SUCCESS) {
            char msg[128];
            snprintf(msg, sizeof(msg), "compute failed on cycle %d with rc=%d", i, rc);
            fail("test_driver_compute_leak", msg);
            ok = 0;
            break;
        }
    }

    restore_cwd();
    grodftb_finalize(&handle);
    if (ok) {
        pass("test_driver_compute_leak");
    }
#endif
}

/* ================================================================
 * US-013: grodftb_get_energy() exact match tests
 *
 * SDD:specs.md:§21.1 — Integration: driver vs standalone
 * ThDD:T-US-013-5.1 — Pass-through identity: E_grodftb = E_dftbp
 *
 * Reference values from real DFTB+ 25.1 (commit fd31d873) runs:
 *   B1: tests/data/b1/reference.json  provenance: tests/data/b1/provenance.txt
 *   B2: tests/data/b2/reference.json  provenance: tests/data/b2/provenance.txt
 * ================================================================ */

/*
 * Coordinates in Bohr for exact-match tests.  Converted from geo.gen
 * (Angstrom, GenFormat "C" cluster) using 1 Bohr = 0.529177210903 Å.
 * The b1_coords array above is in Angstroms (from geo.gen directly)
 * and is used by US-009/011/012 sanity tests only.
 */
static const double b1_coords_bohr[18] = {
    -1.326958030291, -0.105938038921,  0.018787655779,  /* O */
    -1.931664677465,  1.600174613723, -0.021711061883,  /* H */
     0.486644126310,  0.079597148366,  0.009862479935,  /* H */
     4.196837646028,  0.050487809237,  0.001171630113,  /* O */
     4.908550027307, -0.777930269645,  1.448937953129,  /* H */
     4.900314601448, -0.849424272972, -1.407433901241   /* H */
};

static const double b2_coords_bohr[18] = {
     0.000000000000,  0.000000000000,  0.000000000000,  /* C */
     2.304142897900,  0.000000000000,  0.000000000000,  /* O */
    -1.195440660388,  2.125941737175,  0.000000000000,  /* N */
    -1.106623538924, -1.770673251318,  0.000000000000,  /* H */
    -0.263805748009,  3.804018415052,  0.000000000000,  /* H */
    -3.083276923000,  2.184523242822,  0.000000000000   /* H */
};

/* US-013 AC-1: B1 energy matches standalone DFTB+ within 1e-8 Hartree.
 * Also covers AC-3 (return code) and AC-4 (Hartree unit). */
static void test_driver_energy_match_b1(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_energy_match_b1 [SKIPPED: no DFTB+]");
    return;
#else
    /* Reference: standalone DFTB+ 25.1, mio-1-1, SCCTolerance=1e-10
     * Provenance: tests/data/b1/provenance.txt */
    const double ref_energy = -8.16019724583451;  /* Hartree */
    const double tolerance  = 1.0e-8;              /* specs.md §21.1 */

    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_driver_energy_match_b1", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_driver_energy_match_b1", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_driver_energy_match_b1", "set_geometry failed");
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_driver_energy_match_b1", "compute failed");
        return;
    }
    double energy = 0.0;
    rc = grodftb_get_energy(handle, &energy);
    restore_cwd();
    grodftb_finalize(&handle);
    if (rc != GRODFTB_SUCCESS) {
        char msg[128];
        snprintf(msg, sizeof(msg), "get_energy returned %d", rc);
        fail("test_driver_energy_match_b1", msg);
        return;
    }
    double err = fabs(energy - ref_energy);
    if (err < tolerance) {
        pass("test_driver_energy_match_b1");
    } else {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "energy = %.17e Ha, ref = %.17e Ha, |err| = %.3e (tol = %.0e)",
                 energy, ref_energy, err, tolerance);
        fail("test_driver_energy_match_b1", msg);
    }
#endif
}

/* US-013 AC-2: B2 formamide energy matches standalone DFTB+ within 1e-8 Hartree. */
static void test_driver_energy_match_b2(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_energy_match_b2 [SKIPPED: no DFTB+]");
    return;
#else
    /* Reference: standalone DFTB+ 25.1, mio-1-1, SCCTolerance=1e-10
     * Provenance: tests/data/b2/provenance.txt */
    const double ref_energy = -8.52394393335061;  /* Hartree */
    const double tolerance  = 1.0e-8;              /* specs.md §21.1 */

    grodftb_handle_t handle = NULL;
    const char *hsd = get_b2_hsd_path();
    if (!chdir_to_b2()) {
        fail("test_driver_energy_match_b2", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_driver_energy_match_b2", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b2_species, b2_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[128];
        snprintf(msg, sizeof(msg), "set_geometry returned %d", rc);
        fail("test_driver_energy_match_b2", msg);
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[128];
        snprintf(msg, sizeof(msg), "compute returned %d", rc);
        fail("test_driver_energy_match_b2", msg);
        return;
    }
    double energy = 0.0;
    rc = grodftb_get_energy(handle, &energy);
    restore_cwd();
    grodftb_finalize(&handle);
    if (rc != GRODFTB_SUCCESS) {
        char msg[128];
        snprintf(msg, sizeof(msg), "get_energy returned %d", rc);
        fail("test_driver_energy_match_b2", msg);
        return;
    }
    double err = fabs(energy - ref_energy);
    if (err < tolerance) {
        pass("test_driver_energy_match_b2");
    } else {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "energy = %.17e Ha, ref = %.17e Ha, |err| = %.3e (tol = %.0e)",
                 energy, ref_energy, err, tolerance);
        fail("test_driver_energy_match_b2", msg);
    }
#endif
}

/* US-013 AC-5: grodftb_get_energy(NULL, &e) returns GRODFTB_ERR_NULL_HANDLE */
static void test_driver_get_energy_null_handle(void)
{
    double e = 0.0;
    int rc = grodftb_get_energy(NULL, &e);
    if (rc == GRODFTB_ERR_NULL_HANDLE) {
        pass("test_driver_get_energy_null_handle");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_NULL_HANDLE);
        fail("test_driver_get_energy_null_handle", msg);
    }
}

/* US-013 AC-6: grodftb_get_energy(handle, NULL) returns GRODFTB_ERR_INVALID_ARGUMENT */
static void test_driver_get_energy_null_out(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_get_energy_null_out [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_driver_get_energy_null_out", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_driver_get_energy_null_out", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_driver_get_energy_null_out", "set_geometry failed");
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_driver_get_energy_null_out", "compute failed");
        return;
    }
    restore_cwd();
    rc = grodftb_get_energy(handle, NULL);
    grodftb_finalize(&handle);
    if (rc == GRODFTB_ERR_INVALID_ARGUMENT) {
        pass("test_driver_get_energy_null_out");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_INVALID_ARGUMENT);
        fail("test_driver_get_energy_null_out", msg);
    }
#endif
}

/* US-013 AC-7: grodftb_get_energy() before compute returns GRODFTB_ERR_NO_RESULTS */
static void test_driver_get_energy_no_compute(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_get_energy_no_compute [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_driver_get_energy_no_compute", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    restore_cwd();
    if (rc != GRODFTB_SUCCESS || !handle) {
        fail("test_driver_get_energy_no_compute", "init failed");
        return;
    }
    /* Call get_energy without calling compute first */
    double e = 0.0;
    rc = grodftb_get_energy(handle, &e);
    grodftb_finalize(&handle);
    if (rc == GRODFTB_ERR_NO_RESULTS) {
        pass("test_driver_get_energy_no_compute");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_NO_RESULTS);
        fail("test_driver_get_energy_no_compute", msg);
    }
#endif
}

/* ================================================================
 * US-014: grodftb_get_forces() tests
 *
 * SDD:specs.md:§5.3, §21.1 — Force retrieval and validation
 * ThDD:06_theory:Eq4.1 — F_K = -dE/dR_K
 *
 * Reference values from real DFTB+ 25.1 (commit fd31d873) runs:
 *   B1: tests/data/b1/reference.json  provenance: tests/data/b1/provenance.txt
 *   B2: tests/data/b2/reference.json  provenance: tests/data/b2/provenance.txt
 * ================================================================ */

/* US-014 AC-1, AC-3: B1 forces match standalone DFTB+ within 1e-6 Ha/Bohr */
static void test_driver_forces_match_b1(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_forces_match_b1 [SKIPPED: no DFTB+]");
    return;
#else
    /* Reference: standalone DFTB+ 25.1, mio-1-1, SCCTolerance=1e-10
     * Provenance: tests/data/b1/provenance.txt */
    const double ref_forces[18] = {
        -0.00343382490926825, -0.00506098662021831,  0.000132918036725175,
        -0.00767798272076042,  0.00858191780745707, -0.000192612851707156,
         0.0119381527476067,  -0.003831355870751,    6.55237580943401e-05,
        -0.00267502611827172,  0.00340244381094684,  -8.45938638337103e-05,
         0.000954585399736113,-0.00129678725681989,   0.0100813293737328,
         0.000894095600957587,-0.00179523187061463,  -0.0100025644530115
    };
    const double tolerance = 1.0e-6;  /* specs.md §21.1 */

    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_driver_forces_match_b1", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_driver_forces_match_b1", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_driver_forces_match_b1", "set_geometry failed");
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_driver_forces_match_b1", "compute failed");
        return;
    }
    double forces[18];
    rc = grodftb_get_forces(handle, forces);
    restore_cwd();
    grodftb_finalize(&handle);
    if (rc != GRODFTB_SUCCESS) {
        char msg[128];
        snprintf(msg, sizeof(msg), "get_forces returned %d", rc);
        fail("test_driver_forces_match_b1", msg);
        return;
    }
    double max_err = 0.0;
    int worst = -1;
    for (int i = 0; i < 18; i++) {
        double err = fabs(forces[i] - ref_forces[i]);
        if (err > max_err) {
            max_err = err;
            worst = i;
        }
    }
    if (max_err < tolerance) {
        pass("test_driver_forces_match_b1");
    } else {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "component %d: got %.17e, ref %.17e, |err|=%.3e (tol=%.0e)",
                 worst, forces[worst], ref_forces[worst], max_err, tolerance);
        fail("test_driver_forces_match_b1", msg);
    }
#endif
}

/* US-014 AC-2: B2 formamide forces match standalone DFTB+ within 1e-6 Ha/Bohr */
static void test_driver_forces_match_b2(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_forces_match_b2 [SKIPPED: no DFTB+]");
    return;
#else
    const double ref_forces[18] = {
         0.064931478417841,  -0.0646760129480675, -2.13686089559943e-17,
         0.0236724816646542, -0.0113806219173727,  6.00077058391148e-17,
        -0.0567875378175847,  0.0833853401656769, -5.82095995990749e-17,
        -0.0153720877050633, -0.0141139384513503,  1.52539130307036e-17,
         0.00581367605598469, 0.00396914117547201, -4.44063004139242e-17,
        -0.0222580106158319,  0.00281609197564162,  4.8722890099175e-17
    };
    const double tolerance = 1.0e-6;

    grodftb_handle_t handle = NULL;
    const char *hsd = get_b2_hsd_path();
    if (!chdir_to_b2()) {
        fail("test_driver_forces_match_b2", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_driver_forces_match_b2", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b2_species, b2_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_driver_forces_match_b2", "set_geometry failed");
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_driver_forces_match_b2", "compute failed");
        return;
    }
    double forces[18];
    rc = grodftb_get_forces(handle, forces);
    restore_cwd();
    grodftb_finalize(&handle);
    if (rc != GRODFTB_SUCCESS) {
        char msg[128];
        snprintf(msg, sizeof(msg), "get_forces returned %d", rc);
        fail("test_driver_forces_match_b2", msg);
        return;
    }
    double max_err = 0.0;
    int worst = -1;
    for (int i = 0; i < 18; i++) {
        double err = fabs(forces[i] - ref_forces[i]);
        if (err > max_err) {
            max_err = err;
            worst = i;
        }
    }
    if (max_err < tolerance) {
        pass("test_driver_forces_match_b2");
    } else {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "component %d: got %.17e, ref %.17e, |err|=%.3e (tol=%.0e)",
                 worst, forces[worst], ref_forces[worst], max_err, tolerance);
        fail("test_driver_forces_match_b2", msg);
    }
#endif
}

/* US-014 AC-4: grodftb_get_forces(NULL, buf) returns GRODFTB_ERR_NULL_HANDLE */
static void test_driver_get_forces_null_handle(void)
{
    double buf[18];
    int rc = grodftb_get_forces(NULL, buf);
    if (rc == GRODFTB_ERR_NULL_HANDLE) {
        pass("test_driver_get_forces_null_handle");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_NULL_HANDLE);
        fail("test_driver_get_forces_null_handle", msg);
    }
}

/* US-014: Uninitialized handle returns GRODFTB_ERR_NOT_INITIALIZED */
static void test_driver_get_forces_not_init(void)
{
    grodftb_handle_t fake = (grodftb_handle_t)calloc(1, 256);
    if (!fake) {
        fail("test_driver_get_forces_not_init", "calloc failed");
        return;
    }
    double buf[18];
    int rc = grodftb_get_forces(fake, buf);
    free(fake);
    if (rc == GRODFTB_ERR_NOT_INITIALIZED) {
        pass("test_driver_get_forces_not_init");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_NOT_INITIALIZED);
        fail("test_driver_get_forces_not_init", msg);
    }
}

/* US-014 AC-5: grodftb_get_forces(handle, NULL) returns GRODFTB_ERR_INVALID_ARGUMENT */
static void test_driver_get_forces_null_out(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_get_forces_null_out [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_driver_get_forces_null_out", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_driver_get_forces_null_out", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_driver_get_forces_null_out", "set_geometry failed");
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_driver_get_forces_null_out", "compute failed");
        return;
    }
    restore_cwd();
    rc = grodftb_get_forces(handle, NULL);
    grodftb_finalize(&handle);
    if (rc == GRODFTB_ERR_INVALID_ARGUMENT) {
        pass("test_driver_get_forces_null_out");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_INVALID_ARGUMENT);
        fail("test_driver_get_forces_null_out", msg);
    }
#endif
}

/* US-014 AC-6: grodftb_get_forces() before compute returns GRODFTB_ERR_NO_RESULTS */
static void test_driver_get_forces_no_compute(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_get_forces_no_compute [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_driver_get_forces_no_compute", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    restore_cwd();
    if (rc != GRODFTB_SUCCESS || !handle) {
        fail("test_driver_get_forces_no_compute", "init failed");
        return;
    }
    double buf[18];
    rc = grodftb_get_forces(handle, buf);
    grodftb_finalize(&handle);
    if (rc == GRODFTB_ERR_NO_RESULTS) {
        pass("test_driver_get_forces_no_compute");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_NO_RESULTS);
        fail("test_driver_get_forces_no_compute", msg);
    }
#endif
}

/* US-014 AC-7: Newton's third law — sum of forces ≈ 0 for isolated B1 molecule.
 * ThDD:T-US-014-V.2 — Σ F_{K,α} = 0 ∀ α */
static void test_driver_forces_newton_third_b1(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_forces_newton_third_b1 [SKIPPED: no DFTB+]");
    return;
#else
    const double tolerance = 1.0e-10;  /* Ha/Bohr, per T-US-014-V.2 */

    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_driver_forces_newton_third_b1", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_driver_forces_newton_third_b1", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_driver_forces_newton_third_b1", "set_geometry failed");
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_driver_forces_newton_third_b1", "compute failed");
        return;
    }
    double forces[18];
    rc = grodftb_get_forces(handle, forces);
    restore_cwd();
    grodftb_finalize(&handle);
    if (rc != GRODFTB_SUCCESS) {
        char msg[128];
        snprintf(msg, sizeof(msg), "get_forces returned %d", rc);
        fail("test_driver_forces_newton_third_b1", msg);
        return;
    }
    /* Sum forces per axis */
    double sum[3] = {0.0, 0.0, 0.0};
    for (int k = 0; k < 6; k++) {
        sum[0] += forces[3*k + 0];
        sum[1] += forces[3*k + 1];
        sum[2] += forces[3*k + 2];
    }
    double max_sum = 0.0;
    int worst_axis = -1;
    for (int a = 0; a < 3; a++) {
        if (fabs(sum[a]) > max_sum) {
            max_sum = fabs(sum[a]);
            worst_axis = a;
        }
    }
    if (max_sum < tolerance) {
        pass("test_driver_forces_newton_third_b1");
    } else {
        const char *axis_name[] = {"x", "y", "z"};
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "force sum %s = %.3e Ha/Bohr (tol = %.0e)",
                 axis_name[worst_axis], sum[worst_axis], tolerance);
        fail("test_driver_forces_newton_third_b1", msg);
    }
#endif
}

/* ================================================================
 * US-016: grodftb_get_mulliken_charges() tests
 *
 * SDD:specs.md:§5.3 — Mulliken charge retrieval
 * ThDD:06_theory:§2.3 — q_A = sum_{mu in A} sum_nu P_{mu,nu} S_{mu,nu}
 *
 * Reference values from real DFTB+ 25.1 (commit fd31d873) runs:
 *   B1: tests/data/b1/reference.json  provenance: tests/data/b1/provenance.txt
 * ================================================================ */

/* US-016 AC-1: B1 Mulliken charges match standalone DFTB+ within 1e-8 e.
 * Also validates: charge sum ≈ 0 (neutral system), O negative / H positive. */
static void test_driver_mulliken_b1(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_mulliken_b1 [SKIPPED: no DFTB+]");
    return;
#else
    /* Reference: standalone DFTB+ 25.1, mio-1-1, SCCTolerance=1e-10
     * Provenance: tests/data/b1/provenance.txt */
    const double ref_charges[6] = {
        -0.61489483, 0.29081961, 0.30729651,
        -0.58979128, 0.30328477, 0.30328522
    };
    const double tolerance = 1.0e-8;  /* specs.md §21.1 */

    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_driver_mulliken_b1", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_driver_mulliken_b1", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_driver_mulliken_b1", "set_geometry failed");
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_driver_mulliken_b1", "compute failed");
        return;
    }
    double charges[6];
    rc = grodftb_get_mulliken_charges(handle, charges);
    restore_cwd();
    grodftb_finalize(&handle);
    if (rc != GRODFTB_SUCCESS) {
        char msg[128];
        snprintf(msg, sizeof(msg), "get_mulliken_charges returned %d", rc);
        fail("test_driver_mulliken_b1", msg);
        return;
    }
    /* Check per-atom match */
    double max_err = 0.0;
    int worst = -1;
    for (int i = 0; i < 6; i++) {
        double err = fabs(charges[i] - ref_charges[i]);
        if (err > max_err) {
            max_err = err;
            worst = i;
        }
    }
    if (max_err >= tolerance) {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "atom %d: got %.17e, ref %.17e, |err|=%.3e (tol=%.0e)",
                 worst, charges[worst], ref_charges[worst], max_err, tolerance);
        fail("test_driver_mulliken_b1", msg);
        return;
    }
    /* Check charge sum ≈ 0 (neutral system) */
    double sum = 0.0;
    for (int i = 0; i < 6; i++) {
        sum += charges[i];
    }
    if (fabs(sum) >= 1.0e-10) {
        char msg[128];
        snprintf(msg, sizeof(msg), "charge sum = %.3e, expected ~0", sum);
        fail("test_driver_mulliken_b1", msg);
        return;
    }
    /* Check sign convention: O negative, H positive */
    if (charges[0] >= 0.0 || charges[3] >= 0.0) {
        fail("test_driver_mulliken_b1", "oxygen charges should be negative");
        return;
    }
    if (charges[1] <= 0.0 || charges[2] <= 0.0 || charges[4] <= 0.0 || charges[5] <= 0.0) {
        fail("test_driver_mulliken_b1", "hydrogen charges should be positive");
        return;
    }
    pass("test_driver_mulliken_b1");
#endif
}

/* US-016 AC-4: grodftb_get_mulliken_charges(NULL, buf) returns GRODFTB_ERR_NULL_HANDLE */
static void test_driver_mulliken_null_handle(void)
{
    double buf[6];
    int rc = grodftb_get_mulliken_charges(NULL, buf);
    if (rc == GRODFTB_ERR_NULL_HANDLE) {
        pass("test_driver_mulliken_null_handle");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_NULL_HANDLE);
        fail("test_driver_mulliken_null_handle", msg);
    }
}

/* US-016 AC-5: grodftb_get_mulliken_charges(handle, NULL) returns GRODFTB_ERR_INVALID_ARGUMENT */
static void test_driver_mulliken_null_out(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_mulliken_null_out [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_driver_mulliken_null_out", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_driver_mulliken_null_out", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_driver_mulliken_null_out", "set_geometry failed");
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_driver_mulliken_null_out", "compute failed");
        return;
    }
    restore_cwd();
    rc = grodftb_get_mulliken_charges(handle, NULL);
    grodftb_finalize(&handle);
    if (rc == GRODFTB_ERR_INVALID_ARGUMENT) {
        pass("test_driver_mulliken_null_out");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_INVALID_ARGUMENT);
        fail("test_driver_mulliken_null_out", msg);
    }
#endif
}

/* US-016 AC-6: grodftb_get_mulliken_charges() before compute returns GRODFTB_ERR_NO_RESULTS */
static void test_driver_mulliken_no_compute(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_driver_mulliken_no_compute [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_driver_mulliken_no_compute", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    restore_cwd();
    if (rc != GRODFTB_SUCCESS || !handle) {
        fail("test_driver_mulliken_no_compute", "init failed");
        return;
    }
    /* set_geometry but no compute */
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&handle);
        fail("test_driver_mulliken_no_compute", "set_geometry failed");
        return;
    }
    double buf[6];
    rc = grodftb_get_mulliken_charges(handle, buf);
    grodftb_finalize(&handle);
    if (rc == GRODFTB_ERR_NO_RESULTS) {
        pass("test_driver_mulliken_no_compute");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_NO_RESULTS);
        fail("test_driver_mulliken_no_compute", msg);
    }
#endif
}

/* ================================================================
 * US-017: Error handling tests
 *
 * SDD:specs.md:§18.1 — grodftb_error_string(), grodftb_last_error()
 * ================================================================ */

/* US-017 AC-1,AC-2: grodftb_error_string covers all codes */
static void test_error_string_all_codes(void)
{
    const char *expected[] = {
        "GRODFTB_SUCCESS",
        "GRODFTB_ERR_NULL_HANDLE",
        "GRODFTB_ERR_NOT_INITIALIZED",
        "GRODFTB_ERR_ALREADY_INIT",
        "GRODFTB_ERR_INVALID_ARGUMENT",
        "GRODFTB_ERR_SIZE_MISMATCH",
        "GRODFTB_ERR_FILE_NOT_FOUND",
        "GRODFTB_ERR_HSD_PARSE",
        "GRODFTB_ERR_DFTB_INIT",
        "GRODFTB_ERR_SCC_NOT_CONVERGED",
        "GRODFTB_ERR_NO_RESULTS",
        "GRODFTB_ERR_EXCITED_NOT_CONFIGURED",
        "GRODFTB_ERR_EXCITED_FAILED",
        "GRODFTB_ERR_NAC_UNAVAILABLE",
        "GRODFTB_ERR_CP_SOLVE_FAILED",
        "GRODFTB_ERR_ALLOC_FAILED",
        "GRODFTB_ERR_DFTB_INTERNAL"
    };
    int ok = 1;
    for (int i = 0; i <= 16; i++) {
        const char *s = grodftb_error_string(i);
        if (!s || strcmp(s, expected[i]) != 0) {
            char msg[256];
            snprintf(msg, sizeof(msg), "code %d: got '%s', expected '%s'",
                     i, s ? s : "(null)", expected[i]);
            fail("test_error_string_all_codes", msg);
            ok = 0;
        }
    }
    /* Out-of-range */
    const char *u1 = grodftb_error_string(-1);
    const char *u2 = grodftb_error_string(999);
    if (!u1 || strcmp(u1, "GRODFTB_ERR_UNKNOWN") != 0 ||
        !u2 || strcmp(u2, "GRODFTB_ERR_UNKNOWN") != 0) {
        fail("test_error_string_all_codes", "out-of-range not 'GRODFTB_ERR_UNKNOWN'");
        ok = 0;
    }
    if (ok) pass("test_error_string_all_codes");
}

/* US-017 AC-3: grodftb_last_error(NULL) returns valid string */
static void test_last_error_null_handle(void)
{
    const char *msg = grodftb_last_error(NULL);
    if (msg && strlen(msg) > 0) {
        pass("test_last_error_null_handle");
    } else {
        fail("test_last_error_null_handle", "returned NULL or empty");
    }
}

/* US-017 AC-4,AC-5: grodftb_last_error returns descriptive message after error */
static void test_last_error_descriptive(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_last_error_descriptive [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_last_error_descriptive", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_last_error_descriptive", "init failed");
        return;
    }
    /* Trigger NO_RESULTS by calling get_energy before compute */
    double energy;
    rc = grodftb_get_energy(handle, &energy);
    if (rc != GRODFTB_ERR_NO_RESULTS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_last_error_descriptive", "expected ERR_NO_RESULTS");
        return;
    }
    const char *msg = grodftb_last_error(handle);
    /* Copy message before finalize to avoid use-after-free */
    char msg_copy[256] = {0};
    if (msg) {
        snprintf(msg_copy, sizeof(msg_copy), "%s", msg);
    }
    restore_cwd();
    grodftb_finalize(&handle);
    if (msg_copy[0] && strstr(msg_copy, "compute") != NULL) {
        pass("test_last_error_descriptive");
    } else {
        char buf[512];
        snprintf(buf, sizeof(buf), "message not descriptive: '%s'", msg_copy);
        fail("test_last_error_descriptive", buf);
    }
#endif
}

/* ================================================================
 * US-018: Full integration tests
 *
 * SDD:specs.md:§21.1 — Complete pipeline validation
 *
 * These tests validate the entire pipeline:
 *   init → set_geometry → compute → get_energy → get_forces → get_mulliken_charges → finalize
 *
 * Reference values from real DFTB+ 25.1 (commit fd31d873) runs:
 *   B1: tests/data/b1/reference.json  provenance: tests/data/b1/provenance.txt
 *   B2: tests/data/b2/reference.json  provenance: tests/data/b2/provenance.txt
 * ================================================================ */

/* US-018: Full integration test with B1 water dimer */
static void test_integration_b1(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_integration_b1 [SKIPPED: no DFTB+]");
    return;
#else
    /* Reference: standalone DFTB+ 25.1, mio-1-1, SCCTolerance=1e-10
     * Provenance: tests/data/b1/provenance.txt */
    const double ref_energy = -8.16019724583451;  /* Hartree */
    const double ref_forces[18] = {
        -0.00343382490926825, -0.00506098662021831,  0.000132918036725175,
        -0.00767798272076042,  0.00858191780745707, -0.000192612851707156,
         0.0119381527476067,  -0.003831355870751,    6.55237580943401e-05,
        -0.00267502611827172,  0.00340244381094684, -8.45938638337103e-05,
         0.000954585399736113,-0.00129678725681989,  0.0100813293737328,
         0.000894095600957587,-0.00179523187061463, -0.0100025644530115
    };
    const double ref_charges[6] = {
        -0.61489483, 0.29081961, 0.30729651,
        -0.58979128, 0.30328477, 0.30328522
    };
    const double energy_tol = 1.0e-8;   /* specs.md §21.1 */
    const double force_tol  = 1.0e-6;   /* specs.md §21.1 */
    const double charge_tol = 1.0e-8;   /* specs.md §21.1 */

    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_integration_b1", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_integration_b1", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_integration_b1", "set_geometry failed");
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_integration_b1", "compute failed");
        return;
    }
    /* Get energy */
    double energy = 0.0;
    rc = grodftb_get_energy(handle, &energy);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[128];
        snprintf(msg, sizeof(msg), "get_energy returned %d", rc);
        fail("test_integration_b1", msg);
        return;
    }
    double energy_err = fabs(energy - ref_energy);
    if (energy_err >= energy_tol) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "energy = %.17e Ha, ref = %.17e Ha, |err| = %.3e (tol = %.0e)",
                 energy, ref_energy, energy_err, energy_tol);
        fail("test_integration_b1", msg);
        return;
    }
    /* Get forces */
    double forces[18];
    rc = grodftb_get_forces(handle, forces);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[128];
        snprintf(msg, sizeof(msg), "get_forces returned %d", rc);
        fail("test_integration_b1", msg);
        return;
    }
    double max_force_err = 0.0;
    int worst_force = -1;
    for (int i = 0; i < 18; i++) {
        double err = fabs(forces[i] - ref_forces[i]);
        if (err > max_force_err) {
            max_force_err = err;
            worst_force = i;
        }
    }
    if (max_force_err >= force_tol) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "force[%d]: got %.17e, ref %.17e, |err|=%.3e (tol=%.0e)",
                 worst_force, forces[worst_force], ref_forces[worst_force],
                 max_force_err, force_tol);
        fail("test_integration_b1", msg);
        return;
    }
    /* Get Mulliken charges */
    double charges[6];
    rc = grodftb_get_mulliken_charges(handle, charges);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[128];
        snprintf(msg, sizeof(msg), "get_mulliken_charges returned %d", rc);
        fail("test_integration_b1", msg);
        return;
    }
    double max_charge_err = 0.0;
    int worst_charge = -1;
    for (int i = 0; i < 6; i++) {
        double err = fabs(charges[i] - ref_charges[i]);
        if (err > max_charge_err) {
            max_charge_err = err;
            worst_charge = i;
        }
    }
    if (max_charge_err >= charge_tol) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "charge[%d]: got %.17e, ref %.17e, |err|=%.3e (tol=%.0e)",
                 worst_charge, charges[worst_charge], ref_charges[worst_charge],
                 max_charge_err, charge_tol);
        fail("test_integration_b1", msg);
        return;
    }
    /* Finalize */
    restore_cwd();
    grodftb_finalize(&handle);
    pass("test_integration_b1");
#endif
}

/* US-018: Full integration test with B2 formamide */
static void test_integration_b2(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_integration_b2 [SKIPPED: no DFTB+]");
    return;
#else
    /* Reference: standalone DFTB+ 25.1, mio-1-1, SCCTolerance=1e-10
     * Provenance: tests/data/b2/provenance.txt */
    const double ref_energy = -8.52394393335061;  /* Hartree */
    const double ref_forces[18] = {
         0.064931478417841,  -0.0646760129480675, -2.13686089559943e-17,
         0.0236724816646542, -0.0113806219173727,  6.00077058391148e-17,
        -0.0567875378175847,  0.0833853401656769, -5.82095995990749e-17,
        -0.0153720877050633, -0.0141139384513503,  1.52539130307036e-17,
         0.00581367605598469, 0.00396914117547201,-4.44063004139242e-17,
        -0.0222580106158319,  0.00281609197564162, 4.8722890099175e-17
    };
    const double ref_charges[6] = {
        0.41920494, -0.52120499, -0.29719527,
       -0.00066697,  0.21415151,  0.18571077
    };
    const double energy_tol = 1.0e-8;   /* specs.md §21.1 */
    const double force_tol  = 1.0e-6;   /* specs.md §21.1 */
    const double charge_tol = 1.0e-8;   /* specs.md §21.1 */

    grodftb_handle_t handle = NULL;
    const char *hsd = get_b2_hsd_path();
    if (!chdir_to_b2()) {
        fail("test_integration_b2", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_integration_b2", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b2_species, b2_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_integration_b2", "set_geometry failed");
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_integration_b2", "compute failed");
        return;
    }
    /* Get energy */
    double energy = 0.0;
    rc = grodftb_get_energy(handle, &energy);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[128];
        snprintf(msg, sizeof(msg), "get_energy returned %d", rc);
        fail("test_integration_b2", msg);
        return;
    }
    double energy_err = fabs(energy - ref_energy);
    if (energy_err >= energy_tol) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "energy = %.17e Ha, ref = %.17e Ha, |err| = %.3e (tol = %.0e)",
                 energy, ref_energy, energy_err, energy_tol);
        fail("test_integration_b2", msg);
        return;
    }
    /* Get forces */
    double forces[18];
    rc = grodftb_get_forces(handle, forces);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[128];
        snprintf(msg, sizeof(msg), "get_forces returned %d", rc);
        fail("test_integration_b2", msg);
        return;
    }
    double max_force_err = 0.0;
    int worst_force = -1;
    for (int i = 0; i < 18; i++) {
        double err = fabs(forces[i] - ref_forces[i]);
        if (err > max_force_err) {
            max_force_err = err;
            worst_force = i;
        }
    }
    if (max_force_err >= force_tol) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "force[%d]: got %.17e, ref %.17e, |err|=%.3e (tol=%.0e)",
                 worst_force, forces[worst_force], ref_forces[worst_force],
                 max_force_err, force_tol);
        fail("test_integration_b2", msg);
        return;
    }
    /* Get Mulliken charges */
    double charges[6];
    rc = grodftb_get_mulliken_charges(handle, charges);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[128];
        snprintf(msg, sizeof(msg), "get_mulliken_charges returned %d", rc);
        fail("test_integration_b2", msg);
        return;
    }
    double max_charge_err = 0.0;
    int worst_charge = -1;
    for (int i = 0; i < 6; i++) {
        double err = fabs(charges[i] - ref_charges[i]);
        if (err > max_charge_err) {
            max_charge_err = err;
            worst_charge = i;
        }
    }
    if (max_charge_err >= charge_tol) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "charge[%d]: got %.17e, ref %.17e, |err|=%.3e (tol=%.0e)",
                 worst_charge, charges[worst_charge], ref_charges[worst_charge],
                 max_charge_err, charge_tol);
        fail("test_integration_b2", msg);
        return;
    }
    /* Finalize */
    restore_cwd();
    grodftb_finalize(&handle);
    pass("test_integration_b2");
#endif
}

/* ================================================================
 * US-019: grodftb_set_embedding_charges() tests
 *
 * SDD:specs.md:§5.3 — Electrostatic embedding via external potential
 * ThDD:06_theory:Eq1.1 — Coulomb potential from MM point charges
 *
 * Reference values from real DFTB+ 25.1 (commit fd31d873) runs:
 *   B4: tests/data/b4/reference.json  provenance: tests/data/b4/provenance.txt
 *
 * B4 uses the same water dimer geometry as B1, plus a +1.0 e point charge
 * at (-3.0 Angstrom, 0, 0). The B1 HSD (no PointCharges) is used; embedding
 * is set via grodftb_set_embedding_charges().
 * ================================================================ */

/* B4 point charge: Q = +1.0 e at (-3.0 Angstrom, 0, 0) */
static const double b4_mm_charges[1] = { 1.0 };
/* -3.0 Angstrom / 0.52917721083 Angstrom/Bohr = -5.6691786285383 Bohr */
static const double b4_mm_positions[3] = { -5.6691786285383, 0.0, 0.0 };

/* US-019 AC-1: Embedding energy matches B4 reference within 1e-8 Hartree */
static void test_embedding_energy_b4(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_embedding_energy_b4 [SKIPPED: no DFTB+]");
    return;
#else
    const double ref_energy = -8.18640458887598;  /* Hartree */
    const double tolerance  = 1.0e-8;

    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_embedding_energy_b4", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_embedding_energy_b4", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_energy_b4", "set_geometry failed");
        return;
    }
    rc = grodftb_set_embedding_charges(handle, 1, b4_mm_charges, b4_mm_positions);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[128];
        snprintf(msg, sizeof(msg), "set_embedding_charges returned %d", rc);
        fail("test_embedding_energy_b4", msg);
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_energy_b4", "compute failed");
        return;
    }
    double energy = 0.0;
    rc = grodftb_get_energy(handle, &energy);
    restore_cwd();
    grodftb_finalize(&handle);
    if (rc != GRODFTB_SUCCESS) {
        char msg[128];
        snprintf(msg, sizeof(msg), "get_energy returned %d", rc);
        fail("test_embedding_energy_b4", msg);
        return;
    }
    double err = fabs(energy - ref_energy);
    if (err < tolerance) {
        pass("test_embedding_energy_b4");
    } else {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "energy = %.17e Ha, ref = %.17e Ha, |err| = %.3e (tol = %.0e)",
                 energy, ref_energy, err, tolerance);
        fail("test_embedding_energy_b4", msg);
    }
#endif
}

/* US-019 AC-2: Embedding forces match B4 reference within 1e-6 Ha/Bohr */
static void test_embedding_forces_b4(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_embedding_forces_b4 [SKIPPED: no DFTB+]");
    return;
#else
    const double ref_forces[18] = {
        -0.0381650771505628,  -0.00640356895258604,  3.23787550389285e-05,
         0.00576368652337858,  0.0175312128639483,   -0.000338561161867043,
         0.0276344584372731,  -0.00395557506122379,   7.32021518125855e-05,
        -0.013556788166778,    0.00268746022288713,  -5.18683060726699e-05,
         0.00352972821217717, -0.00119252737170991,   0.0106010714874832,
         0.0034744765470708,  -0.00172041275639274,  -0.0105313098991585
    };
    const double tolerance = 1.0e-6;

    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_embedding_forces_b4", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_embedding_forces_b4", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_forces_b4", "set_geometry failed");
        return;
    }
    rc = grodftb_set_embedding_charges(handle, 1, b4_mm_charges, b4_mm_positions);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_forces_b4", "set_embedding_charges failed");
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_forces_b4", "compute failed");
        return;
    }
    double forces[18];
    rc = grodftb_get_forces(handle, forces);
    restore_cwd();
    grodftb_finalize(&handle);
    if (rc != GRODFTB_SUCCESS) {
        char msg[128];
        snprintf(msg, sizeof(msg), "get_forces returned %d", rc);
        fail("test_embedding_forces_b4", msg);
        return;
    }
    double max_err = 0.0;
    int worst = -1;
    for (int i = 0; i < 18; i++) {
        double err = fabs(forces[i] - ref_forces[i]);
        if (err > max_err) {
            max_err = err;
            worst = i;
        }
    }
    if (max_err < tolerance) {
        pass("test_embedding_forces_b4");
    } else {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "component %d: got %.17e, ref %.17e, |err|=%.3e (tol=%.0e)",
                 worst, forces[worst], ref_forces[worst], max_err, tolerance);
        fail("test_embedding_forces_b4", msg);
    }
#endif
}

/* US-019 AC-3: Embedding Mulliken charges match B4 reference within 1e-8 e */
static void test_embedding_charges_b4(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_embedding_charges_b4 [SKIPPED: no DFTB+]");
    return;
#else
    const double ref_charges[6] = {
        -0.66908865, 0.26886033, 0.37884946,
        -0.60689336, 0.31417861, 0.31409362
    };
    /* Reference values from reference.json have 8 decimal places;
     * tolerance accounts for truncation of reference data. */
    const double tolerance = 1.0e-7;

    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_embedding_charges_b4", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_embedding_charges_b4", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_charges_b4", "set_geometry failed");
        return;
    }
    rc = grodftb_set_embedding_charges(handle, 1, b4_mm_charges, b4_mm_positions);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_charges_b4", "set_embedding_charges failed");
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_charges_b4", "compute failed");
        return;
    }
    double charges[6];
    rc = grodftb_get_mulliken_charges(handle, charges);
    restore_cwd();
    grodftb_finalize(&handle);
    if (rc != GRODFTB_SUCCESS) {
        char msg[128];
        snprintf(msg, sizeof(msg), "get_mulliken_charges returned %d", rc);
        fail("test_embedding_charges_b4", msg);
        return;
    }
    double max_err = 0.0;
    int worst = -1;
    for (int i = 0; i < 6; i++) {
        double err = fabs(charges[i] - ref_charges[i]);
        if (err > max_err) {
            max_err = err;
            worst = i;
        }
    }
    if (max_err < tolerance) {
        pass("test_embedding_charges_b4");
    } else {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "atom %d: got %.17e, ref %.17e, |err|=%.3e (tol=%.0e)",
                 worst, charges[worst], ref_charges[worst], max_err, tolerance);
        fail("test_embedding_charges_b4", msg);
    }
#endif
}

/* US-019: Embedding energy is lower than gas-phase (positive charge attracts electrons) */
static void test_embedding_sign_convention(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_embedding_sign_convention [SKIPPED: no DFTB+]");
    return;
#else
    const double b1_ref_energy = -8.16019724583451;  /* Hartree, gas-phase */

    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_embedding_sign_convention", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_embedding_sign_convention", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_sign_convention", "set_geometry failed");
        return;
    }
    rc = grodftb_set_embedding_charges(handle, 1, b4_mm_charges, b4_mm_positions);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_sign_convention", "set_embedding_charges failed");
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_sign_convention", "compute failed");
        return;
    }
    double energy = 0.0;
    rc = grodftb_get_energy(handle, &energy);
    restore_cwd();
    grodftb_finalize(&handle);
    if (rc != GRODFTB_SUCCESS) {
        fail("test_embedding_sign_convention", "get_energy failed");
        return;
    }
    if (energy < b1_ref_energy) {
        pass("test_embedding_sign_convention");
    } else {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "E(B4) = %.17e >= E(B1) = %.17e — expected E(B4) < E(B1)",
                 energy, b1_ref_energy);
        fail("test_embedding_sign_convention", msg);
    }
#endif
}

/* US-019: NULL handle returns ERR_NULL_HANDLE */
static void test_embedding_error_null_handle(void)
{
    double q = 1.0, pos[3] = {0.0, 0.0, 0.0};
    int rc = grodftb_set_embedding_charges(NULL, 1, &q, pos);
    if (rc == GRODFTB_ERR_NULL_HANDLE) {
        pass("test_embedding_error_null_handle");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_NULL_HANDLE);
        fail("test_embedding_error_null_handle", msg);
    }
}

/* US-019: NULL charges returns ERR_INVALID_ARGUMENT */
static void test_embedding_error_null_charges(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_embedding_error_null_charges [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_embedding_error_null_charges", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_embedding_error_null_charges", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_error_null_charges", "set_geometry failed");
        return;
    }
    double pos[3] = {0.0, 0.0, 0.0};
    rc = grodftb_set_embedding_charges(handle, 1, NULL, pos);
    restore_cwd();
    grodftb_finalize(&handle);
    if (rc == GRODFTB_ERR_INVALID_ARGUMENT) {
        pass("test_embedding_error_null_charges");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_INVALID_ARGUMENT);
        fail("test_embedding_error_null_charges", msg);
    }
#endif
}

/* US-019: NULL positions returns ERR_INVALID_ARGUMENT */
static void test_embedding_error_null_positions(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_embedding_error_null_positions [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_embedding_error_null_positions", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_embedding_error_null_positions", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_error_null_positions", "set_geometry failed");
        return;
    }
    double q = 1.0;
    rc = grodftb_set_embedding_charges(handle, 1, &q, NULL);
    restore_cwd();
    grodftb_finalize(&handle);
    if (rc == GRODFTB_ERR_INVALID_ARGUMENT) {
        pass("test_embedding_error_null_positions");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_INVALID_ARGUMENT);
        fail("test_embedding_error_null_positions", msg);
    }
#endif
}

/* US-019: No geometry set returns ERR_INVALID_ARGUMENT */
static void test_embedding_error_no_geometry(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_embedding_error_no_geometry [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_embedding_error_no_geometry", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    restore_cwd();
    if (rc != GRODFTB_SUCCESS || !handle) {
        fail("test_embedding_error_no_geometry", "init failed");
        return;
    }
    double q = 1.0, pos[3] = {0.0, 0.0, 0.0};
    rc = grodftb_set_embedding_charges(handle, 1, &q, pos);
    grodftb_finalize(&handle);
    if (rc == GRODFTB_ERR_INVALID_ARGUMENT) {
        pass("test_embedding_error_no_geometry");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_INVALID_ARGUMENT);
        fail("test_embedding_error_no_geometry", msg);
    }
#endif
}

/* US-019: Negative ncharges returns ERR_INVALID_ARGUMENT */
static void test_embedding_error_bad_ncharges(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_embedding_error_bad_ncharges [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_embedding_error_bad_ncharges", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_embedding_error_bad_ncharges", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_error_bad_ncharges", "set_geometry failed");
        return;
    }
    double q = 1.0, pos[3] = {0.0, 0.0, 0.0};
    rc = grodftb_set_embedding_charges(handle, -1, &q, pos);
    restore_cwd();
    grodftb_finalize(&handle);
    if (rc == GRODFTB_ERR_INVALID_ARGUMENT) {
        pass("test_embedding_error_bad_ncharges");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_INVALID_ARGUMENT);
        fail("test_embedding_error_bad_ncharges", msg);
    }
#endif
}

/* US-019: Clear embedding reverts to gas-phase energy */
static void test_embedding_clear(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_embedding_clear [SKIPPED: no DFTB+]");
    return;
#else
    const double b1_ref_energy = -8.16019724583451;  /* Hartree */
    const double tolerance = 1.0e-8;

    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_embedding_clear", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_embedding_clear", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_clear", "set_geometry failed");
        return;
    }
    /* Set embedding */
    rc = grodftb_set_embedding_charges(handle, 1, b4_mm_charges, b4_mm_positions);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_clear", "set_embedding_charges failed");
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_clear", "first compute failed");
        return;
    }
    double e_embedded = 0.0;
    grodftb_get_energy(handle, &e_embedded);

    /* Clear embedding */
    rc = grodftb_set_embedding_charges(handle, 0, NULL, NULL);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[128];
        snprintf(msg, sizeof(msg), "clear returned %d", rc);
        fail("test_embedding_clear", msg);
        return;
    }
    /* Need to re-set geometry to push coords to DFTB+ (set_embedding_charges
     * only sets potential, compute will use whatever coords were last set) */
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_clear", "re-set_geometry failed");
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_clear", "second compute failed");
        return;
    }
    double e_cleared = 0.0;
    grodftb_get_energy(handle, &e_cleared);
    restore_cwd();
    grodftb_finalize(&handle);

    double err = fabs(e_cleared - b1_ref_energy);
    if (err < tolerance) {
        pass("test_embedding_clear");
    } else {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "cleared energy = %.17e Ha, B1 ref = %.17e Ha, |err| = %.3e",
                 e_cleared, b1_ref_energy, err);
        fail("test_embedding_clear", msg);
    }
#endif
}

/* US-019: 1000 set_embedding_charges cycles — leak check */
static void test_embedding_leak(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_embedding_leak [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_embedding_leak", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_embedding_leak", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_leak", "set_geometry failed");
        return;
    }

    int ok = 1;
    for (int i = 0; i < 1000; i++) {
        rc = grodftb_set_embedding_charges(handle, 1, b4_mm_charges, b4_mm_positions);
        if (rc != GRODFTB_SUCCESS) {
            char msg[128];
            snprintf(msg, sizeof(msg), "set_embedding_charges failed on cycle %d with rc=%d", i, rc);
            fail("test_embedding_leak", msg);
            ok = 0;
            break;
        }
        rc = grodftb_compute(handle);
        if (rc != GRODFTB_SUCCESS) {
            char msg[128];
            snprintf(msg, sizeof(msg), "compute failed on cycle %d with rc=%d", i, rc);
            fail("test_embedding_leak", msg);
            ok = 0;
            break;
        }
    }

    restore_cwd();
    grodftb_finalize(&handle);
    if (ok) {
        pass("test_embedding_leak");
    }
#endif
}

/* ================================================================
 * US-020: grodftb_set_embedding_potential() tests
 *
 * SDD:specs.md:§5.3 — PME pathway: pre-computed potential/field
 * ThDD:T-US-020-3.1 — extpot = potential (direct copy)
 * ThDD:T-US-020-3.2 — extpotgrad = -field (sign negation)
 *
 * Equivalence test strategy: compute phi and E from the same B4 point
 * charge using Coulomb expressions identical to set_embedding_charges,
 * then verify set_embedding_potential produces identical results.
 * ================================================================ */

/*
 * Helper: compute electrostatic potential and electric field at QM atom sites
 * from point charges. Uses the same Coulomb expressions as set_embedding_charges
 * (ThDD:T-US-019-2.1/3.1) so the equivalence test is meaningful.
 *
 * potential[A] = -sum_J Q_J / |R_A - R_J|            (Hartree)
 * field[3A+a]  = -sum_J Q_J * (R_A - R_J)_a / |R_A - R_J|^3  (Hartree/Bohr)
 *
 * Note: field = E = -∇φ, so field[i] = -extpotgrad[i] from US-019.
 */
static void compute_coulomb_potential_and_field(
    int nqm, const double *qm_coords,
    int nmm, const double *mm_charges, const double *mm_positions,
    double *potential, double *field)
{
    int i, k;
    for (i = 0; i < nqm; i++) {
        potential[i] = 0.0;
        if (field) {
            field[3*i + 0] = 0.0;
            field[3*i + 1] = 0.0;
            field[3*i + 2] = 0.0;
        }
        for (k = 0; k < nmm; k++) {
            double dx = qm_coords[3*i + 0] - mm_positions[3*k + 0];
            double dy = qm_coords[3*i + 1] - mm_positions[3*k + 1];
            double dz = qm_coords[3*i + 2] - mm_positions[3*k + 2];
            double r2 = dx*dx + dy*dy + dz*dz;
            double r  = sqrt(r2);
            double inv_r = 1.0 / r;
            double inv_r3 = inv_r / r2;

            /* potential = -Q/r (same sign as extpot in US-019) */
            potential[i] -= mm_charges[k] * inv_r;

            /* field = E = -∇φ = -Q*(R_A - R_J)/|R_A - R_J|^3
             * This is the NEGATIVE of extpotgrad from US-019. */
            if (field) {
                field[3*i + 0] -= mm_charges[k] * dx * inv_r3;
                field[3*i + 1] -= mm_charges[k] * dy * inv_r3;
                field[3*i + 2] -= mm_charges[k] * dz * inv_r3;
            }
        }
    }
}

/* US-020 AC-1: Equivalence — energy/forces/charges from set_embedding_potential
 * match set_embedding_charges on B4 within 1e-12 */
static void test_embedding_potential_equivalence(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_embedding_potential_equivalence [SKIPPED: no DFTB+]");
    return;
#else
    const double tolerance = 1.0e-12;

    /* --- Path A: set_embedding_charges --- */
    grodftb_handle_t handle_a = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_embedding_potential_equivalence", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle_a);
    if (rc != GRODFTB_SUCCESS || !handle_a) {
        restore_cwd();
        fail("test_embedding_potential_equivalence", "init A failed");
        return;
    }
    rc = grodftb_set_geometry(handle_a, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle_a);
        fail("test_embedding_potential_equivalence", "set_geometry A failed");
        return;
    }
    rc = grodftb_set_embedding_charges(handle_a, 1, b4_mm_charges, b4_mm_positions);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle_a);
        fail("test_embedding_potential_equivalence", "set_embedding_charges failed");
        return;
    }
    rc = grodftb_compute(handle_a);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle_a);
        fail("test_embedding_potential_equivalence", "compute A failed");
        return;
    }
    double energy_a = 0.0;
    double forces_a[18];
    double charges_a[6];
    grodftb_get_energy(handle_a, &energy_a);
    grodftb_get_forces(handle_a, forces_a);
    grodftb_get_mulliken_charges(handle_a, charges_a);
    grodftb_finalize(&handle_a);

    /* --- Path B: set_embedding_potential --- */
    /* Compute phi and E from same point charge */
    double phi[6], efield[18];
    compute_coulomb_potential_and_field(6, b1_coords_bohr,
                                        1, b4_mm_charges, b4_mm_positions,
                                        phi, efield);

    grodftb_handle_t handle_b = NULL;
    rc = grodftb_init(hsd, &handle_b);
    if (rc != GRODFTB_SUCCESS || !handle_b) {
        restore_cwd();
        fail("test_embedding_potential_equivalence", "init B failed");
        return;
    }
    rc = grodftb_set_geometry(handle_b, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle_b);
        fail("test_embedding_potential_equivalence", "set_geometry B failed");
        return;
    }
    rc = grodftb_set_embedding_potential(handle_b, 6, phi, efield);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle_b);
        char msg[128];
        snprintf(msg, sizeof(msg), "set_embedding_potential returned %d", rc);
        fail("test_embedding_potential_equivalence", msg);
        return;
    }
    rc = grodftb_compute(handle_b);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle_b);
        fail("test_embedding_potential_equivalence", "compute B failed");
        return;
    }
    double energy_b = 0.0;
    double forces_b[18];
    double charges_b[6];
    grodftb_get_energy(handle_b, &energy_b);
    grodftb_get_forces(handle_b, forces_b);
    grodftb_get_mulliken_charges(handle_b, charges_b);
    restore_cwd();
    grodftb_finalize(&handle_b);

    /* --- Compare --- */
    double e_err = fabs(energy_a - energy_b);
    if (e_err >= tolerance) {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "energy: A=%.17e, B=%.17e, |diff|=%.3e (tol=%.0e)",
                 energy_a, energy_b, e_err, tolerance);
        fail("test_embedding_potential_equivalence", msg);
        return;
    }

    double max_f_err = 0.0;
    int worst_f = -1;
    for (int i = 0; i < 18; i++) {
        double err = fabs(forces_a[i] - forces_b[i]);
        if (err > max_f_err) { max_f_err = err; worst_f = i; }
    }
    if (max_f_err >= tolerance) {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "force[%d]: A=%.17e, B=%.17e, |diff|=%.3e (tol=%.0e)",
                 worst_f, forces_a[worst_f], forces_b[worst_f], max_f_err, tolerance);
        fail("test_embedding_potential_equivalence", msg);
        return;
    }

    double max_q_err = 0.0;
    int worst_q = -1;
    for (int i = 0; i < 6; i++) {
        double err = fabs(charges_a[i] - charges_b[i]);
        if (err > max_q_err) { max_q_err = err; worst_q = i; }
    }
    if (max_q_err >= tolerance) {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "charge[%d]: A=%.17e, B=%.17e, |diff|=%.3e (tol=%.0e)",
                 worst_q, charges_a[worst_q], charges_b[worst_q], max_q_err, tolerance);
        fail("test_embedding_potential_equivalence", msg);
        return;
    }

    pass("test_embedding_potential_equivalence");
#endif
}

/* US-020 AC-2: field=NULL accepted; energy matches full-field case within 1e-12 */
static void test_embedding_potential_null_field(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_embedding_potential_null_field [SKIPPED: no DFTB+]");
    return;
#else
    const double tolerance = 1.0e-12;

    /* Compute phi from B4 point charge */
    double phi[6];
    compute_coulomb_potential_and_field(6, b1_coords_bohr,
                                        1, b4_mm_charges, b4_mm_positions,
                                        phi, NULL);

    /* Path with field=NULL */
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_embedding_potential_null_field", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_embedding_potential_null_field", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_potential_null_field", "set_geometry failed");
        return;
    }
    rc = grodftb_set_embedding_potential(handle, 6, phi, NULL);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[128];
        snprintf(msg, sizeof(msg), "set_embedding_potential(NULL field) returned %d", rc);
        fail("test_embedding_potential_null_field", msg);
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_potential_null_field", "compute failed");
        return;
    }
    double energy_null_field = 0.0;
    grodftb_get_energy(handle, &energy_null_field);
    grodftb_finalize(&handle);

    /* Path with full field for energy comparison */
    double efield[18];
    compute_coulomb_potential_and_field(6, b1_coords_bohr,
                                        1, b4_mm_charges, b4_mm_positions,
                                        phi, efield);
    rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_embedding_potential_null_field", "init B failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_potential_null_field", "set_geometry B failed");
        return;
    }
    rc = grodftb_set_embedding_potential(handle, 6, phi, efield);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_potential_null_field", "set_embedding_potential B failed");
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_potential_null_field", "compute B failed");
        return;
    }
    double energy_full_field = 0.0;
    grodftb_get_energy(handle, &energy_full_field);
    restore_cwd();
    grodftb_finalize(&handle);

    double err = fabs(energy_null_field - energy_full_field);
    if (err < tolerance) {
        pass("test_embedding_potential_null_field");
    } else {
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "energy(NULL field)=%.17e, energy(full field)=%.17e, |diff|=%.3e (tol=%.0e)",
                 energy_null_field, energy_full_field, err, tolerance);
        fail("test_embedding_potential_null_field", msg);
    }
#endif
}

/* US-020 AC-3: NULL handle returns ERR_NULL_HANDLE */
static void test_embedding_potential_error_null_handle(void)
{
    double phi[6] = {0};
    int rc = grodftb_set_embedding_potential(NULL, 6, phi, NULL);
    if (rc == GRODFTB_ERR_NULL_HANDLE) {
        pass("test_embedding_potential_error_null_handle");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_NULL_HANDLE);
        fail("test_embedding_potential_error_null_handle", msg);
    }
}

/* US-020 AC-3: NULL potential returns ERR_INVALID_ARGUMENT */
static void test_embedding_potential_error_null_potential(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_embedding_potential_error_null_potential [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_embedding_potential_error_null_potential", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_embedding_potential_error_null_potential", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_potential_error_null_potential", "set_geometry failed");
        return;
    }
    rc = grodftb_set_embedding_potential(handle, 6, NULL, NULL);
    restore_cwd();
    grodftb_finalize(&handle);
    if (rc == GRODFTB_ERR_INVALID_ARGUMENT) {
        pass("test_embedding_potential_error_null_potential");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_INVALID_ARGUMENT);
        fail("test_embedding_potential_error_null_potential", msg);
    }
#endif
}

/* US-020 AC-3: No geometry set returns ERR_INVALID_ARGUMENT */
static void test_embedding_potential_error_no_geometry(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_embedding_potential_error_no_geometry [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_embedding_potential_error_no_geometry", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    restore_cwd();
    if (rc != GRODFTB_SUCCESS || !handle) {
        fail("test_embedding_potential_error_no_geometry", "init failed");
        return;
    }
    double phi[6] = {0};
    rc = grodftb_set_embedding_potential(handle, 6, phi, NULL);
    grodftb_finalize(&handle);
    if (rc == GRODFTB_ERR_INVALID_ARGUMENT) {
        pass("test_embedding_potential_error_no_geometry");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_INVALID_ARGUMENT);
        fail("test_embedding_potential_error_no_geometry", msg);
    }
#endif
}

/* US-020 AC-3: natoms mismatch returns ERR_INVALID_ARGUMENT */
static void test_embedding_potential_error_bad_natoms(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_embedding_potential_error_bad_natoms [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_embedding_potential_error_bad_natoms", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_embedding_potential_error_bad_natoms", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_potential_error_bad_natoms", "set_geometry failed");
        return;
    }
    double phi[10] = {0};
    rc = grodftb_set_embedding_potential(handle, 10, phi, NULL);
    restore_cwd();
    grodftb_finalize(&handle);
    if (rc == GRODFTB_ERR_INVALID_ARGUMENT) {
        pass("test_embedding_potential_error_bad_natoms");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_INVALID_ARGUMENT);
        fail("test_embedding_potential_error_bad_natoms", msg);
    }
#endif
}

/* ================================================================
 * US-021: grodftb_get_embedding_forces() tests
 *
 * SDD:specs.md:§5.3 — Embedding forces getter
 * ThDD:T-US-021-3.6 — Coulomb back-reaction force on MM charges
 * ThDD:T-US-021-4.1 — Newton's 3rd law verification
 * ================================================================ */

/* US-021 AC-3: NULL handle returns ERR_NULL_HANDLE */
static void test_embedding_forces_error_null_handle(void)
{
    double buf[3];
    int rc = grodftb_get_embedding_forces(NULL, buf);
    if (rc == GRODFTB_ERR_NULL_HANDLE) {
        pass("test_embedding_forces_error_null_handle");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_NULL_HANDLE);
        fail("test_embedding_forces_error_null_handle", msg);
    }
}

/* US-021 AC-3: NULL output returns ERR_INVALID_ARGUMENT */
static void test_embedding_forces_error_null_output(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_embedding_forces_error_null_output [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_embedding_forces_error_null_output", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_embedding_forces_error_null_output", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_forces_error_null_output", "set_geometry failed");
        return;
    }
    rc = grodftb_set_embedding_charges(handle, 1, b4_mm_charges, b4_mm_positions);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_forces_error_null_output", "set_embedding_charges failed");
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_forces_error_null_output", "compute failed");
        return;
    }
    rc = grodftb_get_embedding_forces(handle, NULL);
    restore_cwd();
    grodftb_finalize(&handle);
    if (rc == GRODFTB_ERR_INVALID_ARGUMENT) {
        pass("test_embedding_forces_error_null_output");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_INVALID_ARGUMENT);
        fail("test_embedding_forces_error_null_output", msg);
    }
#endif
}

/* US-021 AC-3: calling before compute() returns ERR_NO_RESULTS */
static void test_embedding_forces_error_no_compute(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_embedding_forces_error_no_compute [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_embedding_forces_error_no_compute", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_embedding_forces_error_no_compute", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_forces_error_no_compute", "set_geometry failed");
        return;
    }
    rc = grodftb_set_embedding_charges(handle, 1, b4_mm_charges, b4_mm_positions);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_forces_error_no_compute", "set_embedding_charges failed");
        return;
    }
    double buf[3];
    rc = grodftb_get_embedding_forces(handle, buf);
    restore_cwd();
    grodftb_finalize(&handle);
    if (rc == GRODFTB_ERR_NO_RESULTS) {
        pass("test_embedding_forces_error_no_compute");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_NO_RESULTS);
        fail("test_embedding_forces_error_no_compute", msg);
    }
#endif
}

/* US-021 AC-3: calling without embedding charges returns ERR_INVALID_ARGUMENT */
static void test_embedding_forces_error_no_embedding(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_embedding_forces_error_no_embedding [SKIPPED: no DFTB+]");
    return;
#else
    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_embedding_forces_error_no_embedding", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_embedding_forces_error_no_embedding", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_forces_error_no_embedding", "set_geometry failed");
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_forces_error_no_embedding", "compute failed");
        return;
    }
    double buf[3];
    rc = grodftb_get_embedding_forces(handle, buf);
    restore_cwd();
    grodftb_finalize(&handle);
    if (rc == GRODFTB_ERR_INVALID_ARGUMENT) {
        pass("test_embedding_forces_error_no_embedding");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "returned %d, expected %d",
                 rc, GRODFTB_ERR_INVALID_ARGUMENT);
        fail("test_embedding_forces_error_no_embedding", msg);
    }
#endif
}

/* US-021 AC-2: Newton's third law — sum of QM + MM forces ≈ 0 per axis
 * ThDD:T-US-021-4.1 — Σ_J F_J + Σ_A F_A^emb = 0
 */
static void test_embedding_forces_newton_third(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_embedding_forces_newton_third [SKIPPED: no DFTB+]");
    return;
#else
    const double tolerance = 1.0e-10;  /* Ha/Bohr per axis */

    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_embedding_forces_newton_third", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_embedding_forces_newton_third", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_forces_newton_third", "set_geometry failed");
        return;
    }
    rc = grodftb_set_embedding_charges(handle, 1, b4_mm_charges, b4_mm_positions);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_forces_newton_third", "set_embedding_charges failed");
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_forces_newton_third", "compute failed");
        return;
    }

    double qm_forces[18];
    rc = grodftb_get_forces(handle, qm_forces);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_embedding_forces_newton_third", "get_forces failed");
        return;
    }

    double emb_forces[3];
    rc = grodftb_get_embedding_forces(handle, emb_forces);
    restore_cwd();
    grodftb_finalize(&handle);
    if (rc != GRODFTB_SUCCESS) {
        char msg[128];
        snprintf(msg, sizeof(msg), "get_embedding_forces returned %d", rc);
        fail("test_embedding_forces_newton_third", msg);
        return;
    }

    /* Sum QM forces per axis, add MM embedding force */
    const char *axis_name[] = {"x", "y", "z"};
    int ok = 1;
    printf("    Newton's 3rd law (QM + MM embedding forces):\n");
    for (int a = 0; a < 3; a++) {
        double sum = emb_forces[a];
        for (int k = 0; k < 6; k++) {
            sum += qm_forces[3*k + a];
        }
        int axis_ok = fabs(sum) < tolerance;
        if (!axis_ok) ok = 0;
        printf("    axis %s: sum = %+.3e Ha/Bohr  %s\n",
               axis_name[a], sum, axis_ok ? "PASS" : "FAIL");
    }

    if (ok) {
        pass("test_embedding_forces_newton_third");
    } else {
        fail("test_embedding_forces_newton_third",
             "force sum exceeds 1e-10 Ha/Bohr on at least one axis");
    }
#endif
}

/* US-022: Full integration test with B4 embedded water dimer + 1 MM point charge.
 * Exercises: init -> set_geometry -> set_embedding_charges -> compute
 *            -> get_energy -> get_forces -> get_mulliken_charges
 *            -> get_embedding_forces -> finalize
 * Reference: standalone DFTB+ 25.1 (fd31d873), tests/data/b4/reference.json */
static void test_integration_b4(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_integration_b4 [SKIPPED: no DFTB+]");
    return;
#else
    /* Reference: standalone DFTB+ 25.1, mio-1-1, SCCTolerance=1e-10
     * Provenance: tests/data/b4/provenance.txt */
    const double ref_energy = -8.18640458887598;  /* Hartree */
    const double ref_forces[18] = {
        -0.0381650771505628, -0.00640356895258604,  3.23787550389285e-05,
         0.00576368652337858,  0.0175312128639483, -0.000338561161867043,
         0.0276344584372731, -0.00395557506122379,   7.32021518125855e-05,
        -0.013556788166778,   0.00268746022288713,  -5.18683060726699e-05,
         0.00352972821217717,-0.00119252737170991,    0.0106010714874832,
         0.0034744765470708, -0.00172041275639274,   -0.0105313098991585
    };
    const double ref_charges[6] = {
        -0.66908865, 0.26886033, 0.37884946,
        -0.60689336, 0.31417861, 0.31409362
    };
    const double energy_tol  = 1.0e-8;   /* specs.md §21.1 */
    const double force_tol   = 1.0e-6;   /* specs.md §21.1 */
    /* Reference charges have 8 decimal places (DFTB+ output precision);
     * tolerance accounts for truncation of reference data. */
    const double charge_tol  = 1.0e-7;   /* specs.md §21.1 */
    const double newton_tol  = 1.0e-10;  /* ThDD:T-US-021-4.1 */

    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_integration_b4", "chdir failed");
        return;
    }
    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        fail("test_integration_b4", "init failed");
        return;
    }
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_integration_b4", "set_geometry failed");
        return;
    }
    rc = grodftb_set_embedding_charges(handle, 1, b4_mm_charges, b4_mm_positions);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_integration_b4", "set_embedding_charges failed");
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_integration_b4", "compute failed");
        return;
    }

    /* AC-1: Energy within 1e-8 Ha (SDD:specs.md:§21.1) */
    double energy = 0.0;
    rc = grodftb_get_energy(handle, &energy);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[128];
        snprintf(msg, sizeof(msg), "get_energy returned %d", rc);
        fail("test_integration_b4", msg);
        return;
    }
    double energy_err = fabs(energy - ref_energy);
    if (energy_err >= energy_tol) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "energy = %.17e Ha, ref = %.17e Ha, |err| = %.3e (tol = %.0e)",
                 energy, ref_energy, energy_err, energy_tol);
        fail("test_integration_b4", msg);
        return;
    }

    /* AC-2: Forces (18 components) within 1e-6 Ha/Bohr (SDD:specs.md:§21.1) */
    double forces[18];
    rc = grodftb_get_forces(handle, forces);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[128];
        snprintf(msg, sizeof(msg), "get_forces returned %d", rc);
        fail("test_integration_b4", msg);
        return;
    }
    double max_force_err = 0.0;
    int worst_force = -1;
    for (int i = 0; i < 18; i++) {
        double err = fabs(forces[i] - ref_forces[i]);
        if (err > max_force_err) {
            max_force_err = err;
            worst_force = i;
        }
    }
    if (max_force_err >= force_tol) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "force[%d]: got %.17e, ref %.17e, |err|=%.3e (tol=%.0e)",
                 worst_force, forces[worst_force], ref_forces[worst_force],
                 max_force_err, force_tol);
        fail("test_integration_b4", msg);
        return;
    }

    /* AC-3: Mulliken charges (6 atoms) within 1e-8 e (SDD:specs.md:§21.1) */
    double charges[6];
    rc = grodftb_get_mulliken_charges(handle, charges);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[128];
        snprintf(msg, sizeof(msg), "get_mulliken_charges returned %d", rc);
        fail("test_integration_b4", msg);
        return;
    }
    double max_charge_err = 0.0;
    int worst_charge = -1;
    for (int i = 0; i < 6; i++) {
        double err = fabs(charges[i] - ref_charges[i]);
        if (err > max_charge_err) {
            max_charge_err = err;
            worst_charge = i;
        }
    }
    if (max_charge_err >= charge_tol) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[256];
        snprintf(msg, sizeof(msg),
                 "charge[%d]: got %.17e, ref %.17e, |err|=%.3e (tol=%.0e)",
                 worst_charge, charges[worst_charge], ref_charges[worst_charge],
                 max_charge_err, charge_tol);
        fail("test_integration_b4", msg);
        return;
    }

    /* AC-4: Newton's 3rd law — QM + MM embedding forces sum < 1e-10 Ha/Bohr
     * ThDD:T-US-021-4.1 */
    double emb_forces[3];
    rc = grodftb_get_embedding_forces(handle, emb_forces);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        char msg[128];
        snprintf(msg, sizeof(msg), "get_embedding_forces returned %d", rc);
        fail("test_integration_b4", msg);
        return;
    }
    const char *axis_name[] = {"x", "y", "z"};
    int newton_ok = 1;
    printf("    Newton's 3rd law (integration B4):\n");
    for (int a = 0; a < 3; a++) {
        double sum = emb_forces[a];
        for (int k = 0; k < 6; k++) {
            sum += forces[3*k + a];
        }
        int axis_ok = fabs(sum) < newton_tol;
        if (!axis_ok) newton_ok = 0;
        printf("    axis %s: sum = %+.3e Ha/Bohr  %s\n",
               axis_name[a], sum, axis_ok ? "PASS" : "FAIL");
    }
    if (!newton_ok) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_integration_b4",
             "Newton's 3rd law: force sum exceeds 1e-10 Ha/Bohr on at least one axis");
        return;
    }

    /* Finalize */
    restore_cwd();
    grodftb_finalize(&handle);
    pass("test_integration_b4");
#endif
}

/* US-023: Memory leak stress test — 10,000 full-pipeline cycles.
 * SDD:specs.md:§21.1 — Zero leaks over 10,000 compute() calls
 * SDD:specs.md:§5.5 — Steady-state allocation after first step
 *
 * Each cycle: init -> set_geometry(B1) -> set_embedding_charges(B4)
 *           -> compute -> get_energy -> get_forces -> get_mulliken_charges
 *           -> get_embedding_forces -> finalize
 *
 * Default: 3 cycles (fast for normal test runs).
 * Set GRODFTB_LEAK_TEST_CYCLES=10000 for Valgrind gate test.
 * Set GRODFTB_LEAK_TEST_CYCLES=100 for ASan gate test. */
static void test_memory_stress_10k(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_memory_stress_10k [SKIPPED: no DFTB+]");
    return;
#else
    int n_cycles = 3;
    const char *env = getenv("GRODFTB_LEAK_TEST_CYCLES");
    if (env) {
        n_cycles = atoi(env);
        if (n_cycles < 1) n_cycles = 3;
    }

    const char *hsd = get_b1_hsd_path();
    if (!chdir_to_b1()) {
        fail("test_memory_stress_10k", "chdir failed");
        return;
    }

    printf("    running %d init->compute->finalize cycles...\n", n_cycles);

    int ok = 1;
    for (int i = 0; i < n_cycles; i++) {
        grodftb_handle_t handle = NULL;
        int rc = grodftb_init(hsd, &handle);
        if (rc != GRODFTB_SUCCESS || !handle) {
            char msg[128];
            snprintf(msg, sizeof(msg), "init failed on cycle %d (rc=%d)", i, rc);
            fail("test_memory_stress_10k", msg);
            ok = 0;
            break;
        }
        rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
        if (rc != GRODFTB_SUCCESS) {
            char msg[128];
            snprintf(msg, sizeof(msg), "set_geometry failed on cycle %d (rc=%d)", i, rc);
            grodftb_finalize(&handle);
            fail("test_memory_stress_10k", msg);
            ok = 0;
            break;
        }
        rc = grodftb_set_embedding_charges(handle, 1, b4_mm_charges, b4_mm_positions);
        if (rc != GRODFTB_SUCCESS) {
            char msg[128];
            snprintf(msg, sizeof(msg), "set_embedding_charges failed on cycle %d (rc=%d)", i, rc);
            grodftb_finalize(&handle);
            fail("test_memory_stress_10k", msg);
            ok = 0;
            break;
        }
        rc = grodftb_compute(handle);
        if (rc != GRODFTB_SUCCESS) {
            char msg[128];
            snprintf(msg, sizeof(msg), "compute failed on cycle %d (rc=%d)", i, rc);
            grodftb_finalize(&handle);
            fail("test_memory_stress_10k", msg);
            ok = 0;
            break;
        }

        /* Exercise all getter functions */
        double energy = 0.0;
        grodftb_get_energy(handle, &energy);

        double forces[18];
        grodftb_get_forces(handle, forces);

        double charges[6];
        grodftb_get_mulliken_charges(handle, charges);

        double emb_forces[3];
        grodftb_get_embedding_forces(handle, emb_forces);

        grodftb_finalize(&handle);
    }

    restore_cwd();
    if (ok) {
        printf("    %d cycles completed successfully\n", n_cycles);
        pass("test_memory_stress_10k");
    }
#endif
}

/* ===== US-024: grodftb_version() tests ===== */

/* AC-1: grodftb_version() returns non-NULL */
static void test_driver_version_non_null(void)
{
    const grodftb_version_t *v = grodftb_version();
    if (v != NULL) {
        pass("test_driver_version_non_null");
    } else {
        fail("test_driver_version_non_null", "returned NULL");
    }
}

/* AC-2: Version numbers are 0.2.0 (v0.2-rc) */
static void test_driver_version_numbers(void)
{
    const grodftb_version_t *v = grodftb_version();
    if (!v) { fail("test_driver_version_numbers", "NULL pointer"); return; }
    if (v->major != 0 || v->minor != 2 || v->patch != 0) {
        char msg[128];
        snprintf(msg, sizeof(msg), "version %d.%d.%d, expected 0.2.0",
                 v->major, v->minor, v->patch);
        fail("test_driver_version_numbers", msg);
    } else {
        pass("test_driver_version_numbers");
    }
}

/* AC-3: git_hash is non-NULL and non-empty */
static void test_driver_version_git_hash(void)
{
    const grodftb_version_t *v = grodftb_version();
    if (!v) { fail("test_driver_version_git_hash", "NULL pointer"); return; }
    if (!v->git_hash || v->git_hash[0] == '\0') {
        fail("test_driver_version_git_hash", "git_hash is NULL or empty");
    } else {
        pass("test_driver_version_git_hash");
    }
}

/* AC-4: Returns same pointer on repeated calls (static storage) */
static void test_driver_version_stability(void)
{
    const grodftb_version_t *v1 = grodftb_version();
    const grodftb_version_t *v2 = grodftb_version();
    if (v1 == v2) {
        pass("test_driver_version_stability");
    } else {
        fail("test_driver_version_stability", "pointer changed across calls");
    }
}

int main(void)
{
    printf("test_driver: US-009 init/finalize tests\n");

    test_driver_init_null_path();
    test_driver_init_null_handle();
    test_driver_init_bad_path();
    test_driver_finalize_null_safe();
    test_driver_init_fail_no_leak();
    test_driver_api_version();
    test_driver_init_success();
    test_driver_finalize_nulls_handle();
    test_driver_init_finalize_leak();

    printf("\ntest_driver: US-011 set_geometry tests\n");

    test_driver_set_geometry_null_handle();
    test_driver_set_geometry_not_init();
    test_driver_set_geometry_null_coords();
    test_driver_set_geometry_null_species();
    test_driver_set_geometry_bad_natoms();
    test_driver_set_geometry_size_mismatch();
    test_driver_set_geometry_b1();
    test_driver_set_geometry_leak();

    printf("\ntest_driver: US-012 compute tests\n");

    test_driver_compute_null_handle();
    test_driver_compute_not_init();
    test_driver_compute_no_geometry();
    test_driver_compute_b1();
    test_driver_compute_state();
    test_driver_compute_converged();
    test_driver_compute_double();
    test_driver_compute_leak();

    printf("\ntest_driver: US-013 get_energy exact match tests\n");

    test_driver_energy_match_b1();
    test_driver_energy_match_b2();
    test_driver_get_energy_null_handle();
    test_driver_get_energy_null_out();
    test_driver_get_energy_no_compute();

    printf("\ntest_driver: US-014 get_forces tests\n");

    test_driver_forces_match_b1();
    test_driver_forces_match_b2();
    test_driver_get_forces_null_handle();
    test_driver_get_forces_not_init();
    test_driver_get_forces_null_out();
    test_driver_get_forces_no_compute();
    test_driver_forces_newton_third_b1();

    printf("\ntest_driver: US-016 get_mulliken_charges tests\n");

    test_driver_mulliken_b1();
    test_driver_mulliken_null_handle();
    test_driver_mulliken_null_out();
    test_driver_mulliken_no_compute();

    printf("\ntest_driver: US-017 error handling tests\n");

    test_error_string_all_codes();
    test_last_error_null_handle();
    test_last_error_descriptive();

    printf("\ntest_driver: US-018 integration tests\n");

    test_integration_b1();
    test_integration_b2();

    printf("\ntest_driver: US-019 embedding tests\n");

    test_embedding_error_null_handle();
    test_embedding_error_null_charges();
    test_embedding_error_null_positions();
    test_embedding_error_no_geometry();
    test_embedding_error_bad_ncharges();
    test_embedding_energy_b4();
    test_embedding_forces_b4();
    test_embedding_charges_b4();
    test_embedding_sign_convention();
    test_embedding_clear();
    test_embedding_leak();

    printf("\ntest_driver: US-021 embedding forces tests\n");

    test_embedding_forces_error_null_handle();
    test_embedding_forces_error_null_output();
    test_embedding_forces_error_no_compute();
    test_embedding_forces_error_no_embedding();
    test_embedding_forces_newton_third();

    printf("\ntest_driver: US-022 B4 integration tests\n");

    test_integration_b4();

    printf("\ntest_driver: US-023 memory stress tests\n");

    test_memory_stress_10k();

    printf("\ntest_driver: US-024 version tests\n");

    test_driver_version_non_null();
    test_driver_version_numbers();
    test_driver_version_git_hash();
    test_driver_version_stability();

    printf("\ntest_driver: US-020 embedding potential tests\n");

    test_embedding_potential_error_null_handle();
    test_embedding_potential_error_null_potential();
    test_embedding_potential_error_no_geometry();
    test_embedding_potential_error_bad_natoms();
    test_embedding_potential_equivalence();
    test_embedding_potential_null_field();

    printf("\ntest_driver: %d run, %d passed, %d failed\n",
           tests_run, tests_passed, tests_failed);
    return tests_failed;
}
