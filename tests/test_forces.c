/*
 * US-015: Finite-difference force validation test for the B1 water dimer.
 *
 * ThDD:06_theory:Eq4.1 — F_K = -dE/dR_K
 * ThDD:T-US-015-2.4 — Central FD: F_FD = -(E(+d) - E(-d)) / (2*d)
 * ThDD:T-US-015-V.1 — |F_analytic - F_FD| < max(1e-4 * |F_analytic|, 1e-6)
 *
 * Both analytic forces and FD forces are computed at runtime from the same
 * DFTB+ handle — no fabricated reference values.
 *
 * SDD:specs.md:§21.1 — FD agreement: relative error < 1e-4 (ground state)
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "grodftb/driver.h"
#include "grodftb/error.h"

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
    chdir(sg_oldcwd);
}

/* B1 water dimer coordinates in Bohr (from geo.gen via 1 Bohr = 0.529177210903 A) */
static const double b1_coords_bohr[18] = {
    -1.326958030291, -0.105938038921,  0.018787655779,  /* O */
    -1.931664677465,  1.600174613723, -0.021711061883,  /* H */
     0.486644126310,  0.079597148366,  0.009862479935,  /* H */
     4.196837646028,  0.050487809237,  0.001171630113,  /* O */
     4.908550027307, -0.777930269645,  1.448937953129,  /* H */
     4.900314601448, -0.849424272972, -1.407433901241   /* H */
};

/* 0-based species: O=0, H=1 */
static const int b1_species[6] = { 0, 1, 1, 0, 1, 1 };

/*
 * US-015: Finite-difference force validation for B1 water dimer.
 *
 * ThDD:06_theory:Eq4.1 — F_K,a = -dE/dR_{K,a}
 * ThDD:T-US-015-2.4 — Central differences: F_FD = -(E(+d) - E(-d)) / (2d)
 * ThDD:T-US-015-V.1 — Combined tolerance: |F_an - F_FD| < max(1e-4*|F_an|, 1e-6)
 *
 * All 37 SCC solves (1 reference + 18*2 perturbed) use the SAME handle.
 * No fabricated reference values — both sides computed at runtime.
 */
static void test_fd_forces_b1(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_fd_forces_b1 [SKIPPED: no DFTB+]");
    return;
#else
    const int natom = 6;
    const int ndof = 3 * natom;  /* 18 */
    const double delta = 1.0e-4; /* Bohr, central differences */

    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();

    if (!chdir_to_b1()) {
        fail("test_fd_forces_b1", "chdir to b1 data dir failed");
        return;
    }

    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        char msg[128];
        snprintf(msg, sizeof(msg), "init returned %d", rc);
        fail("test_fd_forces_b1", msg);
        return;
    }

    /* --- Reference computation: analytic forces at equilibrium geometry --- */
    rc = grodftb_set_geometry(handle, natom, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_fd_forces_b1", "set_geometry (reference) failed");
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_fd_forces_b1", "compute (reference) failed");
        return;
    }

    double f_analytic[18];
    rc = grodftb_get_forces(handle, f_analytic);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_fd_forces_b1", "get_forces (reference) failed");
        return;
    }

    /* --- FD loop: 18 DOFs x 2 perturbations = 36 SCC solves --- */
    /*
     * ThDD:06_theory:Eq4.1 — F_{K,a} = -dE/dR_{K,a}
     * ThDD:T-US-015-2.4 — F_FD = -(E_plus - E_minus) / (2 * delta)
     */
    double f_fd[18];
    int fd_ok = 1;

    for (int i = 0; i < ndof; i++) {
        double coords_pert[18];
        double e_plus = 0.0, e_minus = 0.0;

        /* +delta perturbation */
        memcpy(coords_pert, b1_coords_bohr, sizeof(coords_pert));
        coords_pert[i] += delta;

        rc = grodftb_set_geometry(handle, natom, b1_species, coords_pert);
        if (rc != GRODFTB_SUCCESS) {
            restore_cwd();
            grodftb_finalize(&handle);
            char msg[128];
            snprintf(msg, sizeof(msg), "set_geometry (+d) failed at DOF %d, rc=%d", i, rc);
            fail("test_fd_forces_b1", msg);
            return;
        }
        rc = grodftb_compute(handle);
        if (rc != GRODFTB_SUCCESS) {
            restore_cwd();
            grodftb_finalize(&handle);
            char msg[128];
            snprintf(msg, sizeof(msg), "compute (+d) failed at DOF %d, rc=%d", i, rc);
            fail("test_fd_forces_b1", msg);
            return;
        }
        rc = grodftb_get_energy(handle, &e_plus);
        if (rc != GRODFTB_SUCCESS) {
            restore_cwd();
            grodftb_finalize(&handle);
            char msg[128];
            snprintf(msg, sizeof(msg), "get_energy (+d) failed at DOF %d, rc=%d", i, rc);
            fail("test_fd_forces_b1", msg);
            return;
        }

        /* -delta perturbation */
        memcpy(coords_pert, b1_coords_bohr, sizeof(coords_pert));
        coords_pert[i] -= delta;

        rc = grodftb_set_geometry(handle, natom, b1_species, coords_pert);
        if (rc != GRODFTB_SUCCESS) {
            restore_cwd();
            grodftb_finalize(&handle);
            char msg[128];
            snprintf(msg, sizeof(msg), "set_geometry (-d) failed at DOF %d, rc=%d", i, rc);
            fail("test_fd_forces_b1", msg);
            return;
        }
        rc = grodftb_compute(handle);
        if (rc != GRODFTB_SUCCESS) {
            restore_cwd();
            grodftb_finalize(&handle);
            char msg[128];
            snprintf(msg, sizeof(msg), "compute (-d) failed at DOF %d, rc=%d", i, rc);
            fail("test_fd_forces_b1", msg);
            return;
        }
        rc = grodftb_get_energy(handle, &e_minus);
        if (rc != GRODFTB_SUCCESS) {
            restore_cwd();
            grodftb_finalize(&handle);
            char msg[128];
            snprintf(msg, sizeof(msg), "get_energy (-d) failed at DOF %d, rc=%d", i, rc);
            fail("test_fd_forces_b1", msg);
            return;
        }

        /* ThDD:T-US-015-2.4 — Central difference */
        f_fd[i] = -(e_plus - e_minus) / (2.0 * delta);
    }

    restore_cwd();
    grodftb_finalize(&handle);

    /* --- Compare analytic vs FD forces --- */
    /*
     * ThDD:T-US-015-V.1 — Combined acceptance criterion:
     *   |F_analytic - F_FD| < max(1e-4 * |F_analytic|, 1e-6)
     */
    const char *axis_name[] = {"x", "y", "z"};
    int n_component_pass = 0;
    int n_component_fail = 0;

    printf("    FD force validation (delta = %.0e Bohr):\n", delta);
    printf("    %4s %4s %14s %14s %10s %10s %s\n",
           "atom", "axis", "analytic", "FD", "diff", "rel", "status");

    for (int i = 0; i < ndof; i++) {
        int atom = i / 3;
        int axis = i % 3;
        double diff = fabs(f_analytic[i] - f_fd[i]);
        double rel = (fabs(f_analytic[i]) > 1.0e-15)
                     ? diff / fabs(f_analytic[i])
                     : 0.0;
        double tol = fabs(f_analytic[i]) * 1.0e-4;
        if (tol < 1.0e-6) tol = 1.0e-6;

        int ok = (diff < tol);
        if (ok) {
            n_component_pass++;
        } else {
            n_component_fail++;
            fd_ok = 0;
        }

        printf("    %4d %4s %14.8e %14.8e %10.3e %10.3e %s\n",
               atom, axis_name[axis], f_analytic[i], f_fd[i],
               diff, rel, ok ? "PASS" : "FAIL");
    }

    printf("    Components: %d/%d passed\n", n_component_pass, ndof);

    if (fd_ok) {
        pass("test_fd_forces_b1");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "%d/%d components failed FD validation",
                 n_component_fail, ndof);
        fail("test_fd_forces_b1", msg);
    }
#endif
}

/* B4 point charge: Q = +1.0 e at (-3.0 Angstrom, 0, 0) */
static const double b4_mm_charges[1] = { 1.0 };
/* -3.0 Angstrom / 0.52917721083 Angstrom/Bohr = -5.6691786285383 Bohr */
static const double b4_mm_positions[3] = { -5.6691786285383, 0.0, 0.0 };

/*
 * US-021: Finite-difference validation of embedding forces on MM charges.
 *
 * ThDD:T-US-021-3.6 — F_{J,α} = -Q_J Σ_A Δq_A (R_J - R_A)_α / |R_J - R_A|³
 * ThDD:T-US-021-FD.1 — Central differences: F_FD = -(E(+δ) - E(-δ)) / (2δ)
 * ThDD:T-US-021-3.3 — Hellmann-Feynman: FD includes charge re-polarization,
 *   but at SCC convergence the analytic formula is exact.
 *
 * Perturbs the MM charge position (3 DOFs), not QM atoms.
 * 7 SCC solves total (1 reference + 6 perturbed).
 * Both sides computed at runtime — no fabricated reference values.
 *
 * SDD:specs.md:§21.1 — FD embedding forces: relative error < 10⁻⁴
 */
static void test_fd_embedding_forces_b4(void)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    pass("test_fd_embedding_forces_b4 [SKIPPED: no DFTB+]");
    return;
#else
    const int ndof = 3;           /* x, y, z of 1 MM charge */
    const double delta = 1.0e-4;  /* Bohr, central differences */

    grodftb_handle_t handle = NULL;
    const char *hsd = get_b1_hsd_path();

    if (!chdir_to_b1()) {
        fail("test_fd_embedding_forces_b4", "chdir to b1 data dir failed");
        return;
    }

    int rc = grodftb_init(hsd, &handle);
    if (rc != GRODFTB_SUCCESS || !handle) {
        restore_cwd();
        char msg[128];
        snprintf(msg, sizeof(msg), "init returned %d", rc);
        fail("test_fd_embedding_forces_b4", msg);
        return;
    }

    /* --- Reference computation: analytic embedding forces --- */
    rc = grodftb_set_geometry(handle, 6, b1_species, b1_coords_bohr);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_fd_embedding_forces_b4", "set_geometry failed");
        return;
    }
    rc = grodftb_set_embedding_charges(handle, 1, b4_mm_charges, b4_mm_positions);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_fd_embedding_forces_b4", "set_embedding_charges failed");
        return;
    }
    rc = grodftb_compute(handle);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_fd_embedding_forces_b4", "compute (reference) failed");
        return;
    }

    double f_analytic[3];
    rc = grodftb_get_embedding_forces(handle, f_analytic);
    if (rc != GRODFTB_SUCCESS) {
        restore_cwd();
        grodftb_finalize(&handle);
        fail("test_fd_embedding_forces_b4", "get_embedding_forces failed");
        return;
    }

    /* --- FD loop: 3 DOFs x 2 perturbations = 6 SCC solves --- */
    /*
     * ThDD:T-US-021-FD.1 — Central difference on MM charge position:
     *   F_FD = -(E(R_J + δ e_α) - E(R_J - δ e_α)) / (2δ)
     */
    double f_fd[3];
    int fd_ok = 1;

    for (int i = 0; i < ndof; i++) {
        double pos_pert[3];
        double e_plus = 0.0, e_minus = 0.0;

        /* +delta perturbation */
        memcpy(pos_pert, b4_mm_positions, sizeof(pos_pert));
        pos_pert[i] += delta;

        rc = grodftb_set_embedding_charges(handle, 1, b4_mm_charges, pos_pert);
        if (rc != GRODFTB_SUCCESS) {
            restore_cwd();
            grodftb_finalize(&handle);
            char msg[128];
            snprintf(msg, sizeof(msg), "set_embedding_charges (+d) failed at DOF %d", i);
            fail("test_fd_embedding_forces_b4", msg);
            return;
        }
        rc = grodftb_compute(handle);
        if (rc != GRODFTB_SUCCESS) {
            restore_cwd();
            grodftb_finalize(&handle);
            char msg[128];
            snprintf(msg, sizeof(msg), "compute (+d) failed at DOF %d", i);
            fail("test_fd_embedding_forces_b4", msg);
            return;
        }
        rc = grodftb_get_energy(handle, &e_plus);
        if (rc != GRODFTB_SUCCESS) {
            restore_cwd();
            grodftb_finalize(&handle);
            char msg[128];
            snprintf(msg, sizeof(msg), "get_energy (+d) failed at DOF %d", i);
            fail("test_fd_embedding_forces_b4", msg);
            return;
        }

        /* -delta perturbation */
        memcpy(pos_pert, b4_mm_positions, sizeof(pos_pert));
        pos_pert[i] -= delta;

        rc = grodftb_set_embedding_charges(handle, 1, b4_mm_charges, pos_pert);
        if (rc != GRODFTB_SUCCESS) {
            restore_cwd();
            grodftb_finalize(&handle);
            char msg[128];
            snprintf(msg, sizeof(msg), "set_embedding_charges (-d) failed at DOF %d", i);
            fail("test_fd_embedding_forces_b4", msg);
            return;
        }
        rc = grodftb_compute(handle);
        if (rc != GRODFTB_SUCCESS) {
            restore_cwd();
            grodftb_finalize(&handle);
            char msg[128];
            snprintf(msg, sizeof(msg), "compute (-d) failed at DOF %d", i);
            fail("test_fd_embedding_forces_b4", msg);
            return;
        }
        rc = grodftb_get_energy(handle, &e_minus);
        if (rc != GRODFTB_SUCCESS) {
            restore_cwd();
            grodftb_finalize(&handle);
            char msg[128];
            snprintf(msg, sizeof(msg), "get_energy (-d) failed at DOF %d", i);
            fail("test_fd_embedding_forces_b4", msg);
            return;
        }

        /* ThDD:T-US-021-FD.1 — Central difference */
        f_fd[i] = -(e_plus - e_minus) / (2.0 * delta);
    }

    restore_cwd();
    grodftb_finalize(&handle);

    /* --- Compare analytic vs FD embedding forces --- */
    /*
     * Tolerance: max(1e-4 * |F_analytic|, 1e-6) Ha/Bohr
     * From specs.md §21.1 and docs/theory/US-021/04_verification_criteria.md
     */
    const char *axis_name[] = {"x", "y", "z"};
    int n_pass = 0, n_fail = 0;

    printf("    FD embedding force validation (delta = %.0e Bohr):\n", delta);
    printf("    %4s %14s %14s %10s %10s %s\n",
           "axis", "analytic", "FD", "diff", "rel", "status");

    for (int i = 0; i < ndof; i++) {
        double diff = fabs(f_analytic[i] - f_fd[i]);
        double rel = (fabs(f_analytic[i]) > 1.0e-15)
                     ? diff / fabs(f_analytic[i])
                     : 0.0;
        double tol = fabs(f_analytic[i]) * 1.0e-4;
        if (tol < 1.0e-6) tol = 1.0e-6;

        int ok = (diff < tol);
        if (ok) {
            n_pass++;
        } else {
            n_fail++;
            fd_ok = 0;
        }

        printf("    %4s %14.8e %14.8e %10.3e %10.3e %s\n",
               axis_name[i], f_analytic[i], f_fd[i],
               diff, rel, ok ? "PASS" : "FAIL");
    }

    printf("    Components: %d/%d passed\n", n_pass, ndof);

    if (fd_ok) {
        pass("test_fd_embedding_forces_b4");
    } else {
        char msg[128];
        snprintf(msg, sizeof(msg), "%d/%d components failed FD validation",
                 n_fail, ndof);
        fail("test_fd_embedding_forces_b4", msg);
    }
#endif
}

int main(void)
{
    printf("test_forces: US-015 finite-difference force validation\n");

    test_fd_forces_b1();

    printf("\ntest_forces: US-021 finite-difference embedding force validation\n");

    test_fd_embedding_forces_b4();

    printf("\ntest_forces: %d run, %d passed, %d failed\n",
           tests_run, tests_passed, tests_failed);
    return tests_failed;
}
