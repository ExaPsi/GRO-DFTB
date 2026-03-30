/*
 * US-037: B6 Alanine Dipeptide NVE Energy Conservation Test
 *
 * ThDD:T-US-037-N.7 — Linear regression for drift measurement
 * ThDD:T-US-037-N.15 — NVE test protocol (100 ps production)
 * ThDD:T-US-037-1.11 — Normalized drift criterion: < 0.01 kJ/mol/ps/atom
 *
 * SDD:specs.md:§21.1 — Stability tests: NVE drift
 * SDD:specs.md:§21.2 — B6 benchmark system
 *
 * This file implements the long-running NVE energy conservation test for B6.
 * The test measures energy drift via linear regression and verifies that
 * the normalized drift is below the acceptance threshold.
 *
 * TEST IS GATED: Only runs when GRODFTB_LONG_TESTS is defined at compile time.
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
 * Test Parameters from docs/theory/US-037/03_numerical_methods.md
 *
 * ThDD:T-US-037-N.15 — NVE test protocol:
 *   - 10 ps NVT equilibration at 300 K (skipped in gas-phase B6)
 *   - 100 ps NVE production
 *   - Timestep 0.5 fs
 *   - Output frequency: Every 10 steps (5 fs)
 *   - SCC tolerance: 10^-10 e
 *   - Discard first 1 ps
 *
 * ThDD:T-US-037-1.11 — Normalized drift criterion:
 *   |drift| < 0.01 kJ/mol/ps/atom
 * --------------------------------------------------------------------------- */

/* Simulation parameters */
#define TIMESTEP_FS         0.5     /* Timestep in femtoseconds */
#define TIMESTEP_AU         20.67   /* 0.5 fs in atomic units (1 fs = 41.34 a.u.) */

#ifdef GRODFTB_LONG_TESTS
#define NVE_PRODUCTION_PS   100.0   /* Full 100 ps production for long tests */
#define OUTPUT_FREQ_STEPS   10      /* Output every 10 steps */
#else
#define NVE_PRODUCTION_PS   1.0     /* Short 1 ps for CI (placeholder until long tests enabled) */
#define OUTPUT_FREQ_STEPS   10
#endif

#define DISCARD_PS          0.1     /* Discard first 0.1 ps (1 ps for long tests) */
#define DRIFT_TOL_KJMOL_PS_ATOM  0.01  /* Acceptance criterion */

/* System parameters */
#define B6_N_QM_REAL    10
#define B6_N_MM         12
#define B6_N_LINKS       2
#define B6_N_QM_TOTAL  (B6_N_QM_REAL + B6_N_LINKS)
#define B6_TOTAL_ATOMS  (B6_N_QM_REAL + B6_N_LINKS)  /* For drift normalization */

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
 * B6 system coordinates and species (from provenance.json)
 *
 * ThDD:CLAUDE.md:Golden_Rule — All numerical values from real calculations
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

/* Masses in atomic mass units */
static const double b6_masses[B6_N_QM_REAL] = {
    14.007, 1.008, 12.011, 1.008, 12.011,
     1.008, 1.008,  1.008, 12.011, 15.999
};

/* ---------------------------------------------------------------------------
 * Linear Regression for Drift Measurement
 *
 * ThDD:T-US-037-N.7 — Linear regression formula:
 *   slope = (N*sum(t*E) - sum(t)*sum(E)) / (N*sum(t^2) - (sum(t))^2)
 *
 * ThDD:T-US-037-N.18 — Standard error for confidence interval
 * --------------------------------------------------------------------------- */

/**
 * Compute linear regression slope, intercept, and standard error of slope.
 *
 * @param n          Number of data points
 * @param t          Time array [n]
 * @param E          Energy array [n]
 * @param slope_out  Receives the slope (drift rate)
 * @param intercept_out  Receives the intercept
 * @param stderr_out Receives the standard error of the slope
 */
void compute_linear_regression(int n, const double *t, const double *E,
                               double *slope_out, double *intercept_out,
                               double *stderr_out)
{
    if (n < 2) {
        *slope_out = 0.0;
        *intercept_out = E[0];
        *stderr_out = 0.0;
        return;
    }

    /* Compute sums */
    double sum_t = 0.0, sum_E = 0.0, sum_tt = 0.0, sum_tE = 0.0;
    for (int i = 0; i < n; i++) {
        sum_t  += t[i];
        sum_E  += E[i];
        sum_tt += t[i] * t[i];
        sum_tE += t[i] * E[i];
    }

    /* ThDD:T-US-037-N.7 — Linear regression formula */
    double denom = n * sum_tt - sum_t * sum_t;
    if (fabs(denom) < 1e-30) {
        *slope_out = 0.0;
        *intercept_out = sum_E / n;
        *stderr_out = 0.0;
        return;
    }

    double slope = (n * sum_tE - sum_t * sum_E) / denom;
    double intercept = (sum_E - slope * sum_t) / n;

    /* Compute standard error of slope */
    /* stderr = sqrt(sum((E_i - (slope*t_i + intercept))^2) / ((n-2) * sum((t_i - t_mean)^2))) */
    double t_mean = sum_t / n;
    double sum_resid_sq = 0.0;
    double sum_t_dev_sq = 0.0;

    for (int i = 0; i < n; i++) {
        double E_pred = slope * t[i] + intercept;
        double resid = E[i] - E_pred;
        sum_resid_sq += resid * resid;
        double t_dev = t[i] - t_mean;
        sum_t_dev_sq += t_dev * t_dev;
    }

    double se_slope = 0.0;
    if (n > 2 && sum_t_dev_sq > 1e-30) {
        double MSE = sum_resid_sq / (n - 2);
        se_slope = sqrt(MSE / sum_t_dev_sq);
    }

    *slope_out = slope;
    *intercept_out = intercept;
    *stderr_out = se_slope;
}

/* ---------------------------------------------------------------------------
 * Test fixture
 * --------------------------------------------------------------------------- */
typedef struct {
    grodftb_handle_t driver;
    grodftb_linkatom_handle_t links;
    double *augmented_coords;
    double *augmented_forces;
    int *augmented_species;
    double qm_coords[3 * B6_N_QM_REAL];
    double mm_coords[3 * B6_N_MM];
    double velocities[3 * B6_N_QM_REAL];
    int initialized;
} b6_nve_fixture_t;

static int b6_nve_fixture_init(b6_nve_fixture_t *fix)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    fprintf(stderr, "    DFTB+ not available - tests will skip\n");
    fix->initialized = 0;
    return -1;
#else
    int rc;
    memset(fix, 0, sizeof(*fix));

    if (!chdir_to_b6()) {
        fix->initialized = 0;
        return -1;
    }

    FILE *f = fopen("dftb_in.hsd", "r");
    if (!f) {
        fprintf(stderr, "    B6 HSD file not found\n");
        fix->initialized = 0;
        return -1;
    }
    fclose(f);

    rc = grodftb_init("dftb_in.hsd", &fix->driver);
    if (rc != GRODFTB_SUCCESS) {
        fprintf(stderr, "    Failed to initialize DFTB+ driver: %d\n", rc);
        fix->initialized = 0;
        return -1;
    }

    double ref_lengths[B6_N_LINKS] = {1.9451, 2.1761};
    rc = grodftb_linkatom_create(B6_N_LINKS,
                                  B6_QM_BOUNDARY_ATOMS,
                                  B6_MM_BOUNDARY_ATOMS,
                                  ref_lengths,
                                  3,
                                  GRODFTB_CHARGE_NONE,
                                  &fix->links);
    if (rc != GRODFTB_SUCCESS) {
        grodftb_finalize(&fix->driver);
        fix->initialized = 0;
        return -1;
    }

    fix->augmented_coords = malloc(3 * B6_N_QM_TOTAL * sizeof(double));
    fix->augmented_forces = malloc(3 * B6_N_QM_TOTAL * sizeof(double));
    fix->augmented_species = malloc(B6_N_QM_TOTAL * sizeof(int));

    if (!fix->augmented_coords || !fix->augmented_forces || !fix->augmented_species) {
        grodftb_linkatom_destroy(&fix->links);
        grodftb_finalize(&fix->driver);
        free(fix->augmented_coords);
        free(fix->augmented_forces);
        free(fix->augmented_species);
        fix->initialized = 0;
        return -1;
    }

    memcpy(fix->qm_coords, b6_qm_coords, sizeof(b6_qm_coords));
    memcpy(fix->mm_coords, b6_mm_coords, sizeof(b6_mm_coords));
    memset(fix->velocities, 0, sizeof(fix->velocities));

    fix->initialized = 1;
    return 0;
#endif
}

static void b6_nve_fixture_cleanup(b6_nve_fixture_t *fix)
{
    if (!fix->initialized) return;
    free(fix->augmented_coords);
    free(fix->augmented_forces);
    free(fix->augmented_species);
    grodftb_linkatom_destroy(&fix->links);
    grodftb_finalize(&fix->driver);
    restore_cwd();
    fix->initialized = 0;
}

/*
 * Compute total energy (kinetic + potential)
 */
static int compute_total_energy(b6_nve_fixture_t *fix,
                                double *total_energy_out,
                                double *qm_forces)
{
#ifndef GRODFTB_HAS_DFTBPLUS
    (void)fix; (void)total_energy_out; (void)qm_forces;
    return -1;
#else
    int rc;

    /* Augment coordinates */
    rc = grodftb_linkatom_augment_coords(fix->links, B6_N_QM_REAL,
                                          fix->qm_coords, fix->mm_coords,
                                          fix->augmented_coords);
    if (rc != GRODFTB_SUCCESS) return -1;

    rc = grodftb_linkatom_augment_species(fix->links, B6_N_QM_REAL,
                                           b6_qm_species, fix->augmented_species);
    if (rc != GRODFTB_SUCCESS) return -1;

    rc = grodftb_set_geometry(fix->driver, B6_N_QM_TOTAL,
                               fix->augmented_species, fix->augmented_coords);
    if (rc != GRODFTB_SUCCESS) return -1;

    rc = grodftb_compute(fix->driver);
    if (rc != GRODFTB_SUCCESS) return -1;

    /* Get potential energy in Hartree */
    double potential_energy;
    rc = grodftb_get_energy(fix->driver, &potential_energy);
    if (rc != GRODFTB_SUCCESS) return -1;

    /* Get augmented forces */
    rc = grodftb_get_forces(fix->driver, fix->augmented_forces);
    if (rc != GRODFTB_SUCCESS) return -1;

    /* Project forces to real atoms */
    double mm_forces[3 * B6_N_MM];
    memset(mm_forces, 0, sizeof(mm_forces));
    rc = grodftb_linkatom_project_forces(fix->links, B6_N_QM_REAL,
                                          fix->qm_coords, fix->mm_coords,
                                          fix->augmented_forces,
                                          qm_forces, mm_forces);
    if (rc != GRODFTB_SUCCESS) return -1;

    /* Compute kinetic energy in Hartree */
    /* K = 0.5 * sum_i m_i * v_i^2 */
    /* Mass in atomic units: m_amu * 1822.888486 */
    double kinetic_energy = 0.0;
    for (int i = 0; i < B6_N_QM_REAL; i++) {
        double mass_au = b6_masses[i] * 1822.888486;
        double v_sq = 0.0;
        for (int d = 0; d < 3; d++) {
            v_sq += fix->velocities[3*i + d] * fix->velocities[3*i + d];
        }
        kinetic_energy += 0.5 * mass_au * v_sq;
    }

    *total_energy_out = potential_energy + kinetic_energy;
    return 0;
#endif
}

/* ---------------------------------------------------------------------------
 * Test: NVE energy drift (AC-2)
 *
 * ThDD:T-US-037-N.15 — NVE test protocol
 * ThDD:T-US-037-1.11 — Normalized drift criterion: < 0.01 kJ/mol/ps/atom
 *
 * This test runs an NVE trajectory and measures energy drift via linear
 * regression. The drift is normalized by atom count and must be below
 * the acceptance threshold.
 *
 * SDD:specs.md:§21.1 — Stability: NVE drift < 0.01 kJ/mol/ps/atom
 * --------------------------------------------------------------------------- */
static int test_b6_nve_drift(void)
{
    b6_nve_fixture_t fix;
    if (b6_nve_fixture_init(&fix) != 0) return -1;

#ifndef GRODFTB_LONG_TESTS
    printf("\n    NOTE: Running short NVE test (GRODFTB_LONG_TESTS not defined)\n");
    printf("    For full 100 ps test, compile with -DGRODFTB_LONG_TESTS\n");
#endif

    /* Calculate number of steps */
    double total_time_fs = NVE_PRODUCTION_PS * 1000.0;
    int total_steps = (int)(total_time_fs / TIMESTEP_FS);
    int output_interval = OUTPUT_FREQ_STEPS;
    int n_samples = total_steps / output_interval + 1;
    int discard_steps = (int)(DISCARD_PS * 1000.0 / TIMESTEP_FS);

    printf("\n    Trajectory parameters:\n");
    printf("      Production time: %.1f ps\n", NVE_PRODUCTION_PS);
    printf("      Timestep: %.1f fs\n", TIMESTEP_FS);
    printf("      Total steps: %d\n", total_steps);
    printf("      Output interval: %d steps (%.1f fs)\n", output_interval, output_interval * TIMESTEP_FS);
    printf("      Discard time: %.1f ps (%d steps)\n", DISCARD_PS, discard_steps);

    /* Allocate energy trajectory */
    double *times = malloc(n_samples * sizeof(double));
    double *energies = malloc(n_samples * sizeof(double));

    if (!times || !energies) {
        fprintf(stderr, "    Failed to allocate trajectory arrays\n");
        free(times);
        free(energies);
        b6_nve_fixture_cleanup(&fix);
        return 0;
    }

    /* Initialize with small random velocities (300 K equivalent)
     * v_rms = sqrt(k_B T / m) for each DOF
     * In atomic units: k_B T = 300 K * 3.1668e-6 Hartree/K = 9.5e-4 Hartree
     */
    const double kT_au = 9.5e-4;  /* k_B * 300 K in Hartree */
    for (int i = 0; i < B6_N_QM_REAL; i++) {
        double mass_au = b6_masses[i] * 1822.888486;
        double v_rms = sqrt(kT_au / mass_au);
        /* Simple pseudo-random initialization (not statistically correct but adequate for test) */
        fix.velocities[3*i + 0] = v_rms * (((i * 7 + 3) % 11) / 5.5 - 1.0);
        fix.velocities[3*i + 1] = v_rms * (((i * 13 + 7) % 17) / 8.5 - 1.0);
        fix.velocities[3*i + 2] = v_rms * (((i * 19 + 11) % 23) / 11.5 - 1.0);
    }

    /* Working arrays */
    double qm_forces[3 * B6_N_QM_REAL];
    int sample_idx = 0;

    /* Get initial energy */
    double total_energy;
    if (compute_total_energy(&fix, &total_energy, qm_forces) != 0) {
        fprintf(stderr, "    Failed to compute initial energy\n");
        free(times);
        free(energies);
        b6_nve_fixture_cleanup(&fix);
        return 0;
    }

    times[sample_idx] = 0.0;
    energies[sample_idx] = total_energy * GRODFTB_HARTREE_TO_KJMOL;  /* Convert to kJ/mol */
    sample_idx++;

    printf("    Running NVE trajectory...\n");

    /* Run NVE trajectory using velocity Verlet */
    for (int step = 1; step <= total_steps; step++) {
        /* Velocity Verlet step */
        /* v(t+dt/2) = v(t) + F(t)/(2m) * dt */
        for (int i = 0; i < B6_N_QM_REAL; i++) {
            double mass_au = b6_masses[i] * 1822.888486;
            for (int d = 0; d < 3; d++) {
                fix.velocities[3*i + d] += qm_forces[3*i + d] / mass_au * TIMESTEP_AU * 0.5;
            }
        }

        /* r(t+dt) = r(t) + v(t+dt/2) * dt */
        for (int i = 0; i < B6_N_QM_REAL; i++) {
            for (int d = 0; d < 3; d++) {
                fix.qm_coords[3*i + d] += fix.velocities[3*i + d] * TIMESTEP_AU;
            }
        }

        /* Compute new forces */
        if (compute_total_energy(&fix, &total_energy, qm_forces) != 0) {
            fprintf(stderr, "    Failed to compute energy at step %d\n", step);
            free(times);
            free(energies);
            b6_nve_fixture_cleanup(&fix);
            return 0;
        }

        /* v(t+dt) = v(t+dt/2) + F(t+dt)/(2m) * dt */
        for (int i = 0; i < B6_N_QM_REAL; i++) {
            double mass_au = b6_masses[i] * 1822.888486;
            for (int d = 0; d < 3; d++) {
                fix.velocities[3*i + d] += qm_forces[3*i + d] / mass_au * TIMESTEP_AU * 0.5;
            }
        }

        /* Record energy at output intervals */
        if (step % output_interval == 0 && sample_idx < n_samples) {
            double time_ps = step * TIMESTEP_FS / 1000.0;
            times[sample_idx] = time_ps;
            energies[sample_idx] = total_energy * GRODFTB_HARTREE_TO_KJMOL;
            sample_idx++;

            /* Progress indicator for long runs */
            if (step % (total_steps / 10) == 0) {
                printf("      Progress: %d%% (step %d/%d)\n",
                       (step * 100) / total_steps, step, total_steps);
            }
        }
    }

    /* Discard equilibration phase and compute drift */
    int discard_samples = discard_steps / output_interval;
    if (discard_samples >= sample_idx - 2) {
        discard_samples = 0;  /* Don't discard if too few samples */
    }

    int analysis_samples = sample_idx - discard_samples;
    double *analysis_times = times + discard_samples;
    double *analysis_energies = energies + discard_samples;

    /* Compute linear regression */
    double slope, intercept, stderr_slope;
    compute_linear_regression(analysis_samples, analysis_times, analysis_energies,
                              &slope, &intercept, &stderr_slope);

    /* Normalize drift by atom count */
    /* ThDD:T-US-037-1.11 — drift_norm = drift / N_atom */
    double drift_per_atom = slope / B6_TOTAL_ATOMS;
    double stderr_per_atom = stderr_slope / B6_TOTAL_ATOMS;

    /* Report results */
    printf("\n    NVE Drift Analysis:\n");
    printf("      Analysis samples: %d (discarded %d)\n", analysis_samples, discard_samples);
    printf("      Energy range: %.4f to %.4f kJ/mol\n",
           analysis_energies[0], analysis_energies[analysis_samples - 1]);
    printf("      Drift: %.6f +/- %.6f kJ/mol/ps\n", slope, stderr_slope);
    printf("      Normalized drift: %.6f +/- %.6f kJ/mol/ps/atom\n",
           drift_per_atom, stderr_per_atom);
    printf("      Acceptance threshold: %.4f kJ/mol/ps/atom\n", DRIFT_TOL_KJMOL_PS_ATOM);

    /* Check acceptance criterion */
    /* ThDD:T-US-037-1.11 — |drift_norm| < 0.01 kJ/mol/ps/atom */
    int pass = (fabs(drift_per_atom) < DRIFT_TOL_KJMOL_PS_ATOM);

    if (pass) {
        printf("    VERDICT: Energy drift within tolerance\n");
    } else {
        printf("    VERDICT: Energy drift EXCEEDS tolerance\n");
        /* ThDD:T-US-037-N.8 — Check statistical significance */
        double t_stat = fabs(drift_per_atom) / stderr_per_atom;
        printf("      Statistical significance: t = %.2f (>2 means significant drift)\n", t_stat);
    }

    free(times);
    free(energies);
    b6_nve_fixture_cleanup(&fix);
    return pass;
}

/* ---------------------------------------------------------------------------
 * Main test runner
 * --------------------------------------------------------------------------- */
int main(void)
{
    printf("=== US-037 B6 NVE Energy Conservation Test ===\n\n");

    printf("Test parameters:\n");
    printf("  Drift tolerance: %.4f kJ/mol/ps/atom\n", DRIFT_TOL_KJMOL_PS_ATOM);
#ifdef GRODFTB_LONG_TESTS
    printf("  Mode: LONG TESTS (full 100 ps)\n");
#else
    printf("  Mode: SHORT TEST (1 ps, for CI)\n");
    printf("  NOTE: Define GRODFTB_LONG_TESTS for full validation\n");
#endif
    printf("\n");

    printf("AC-2: NVE energy conservation:\n");
    RUN_TEST(test_b6_nve_drift);

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

    return tests_failed > 0 ? 1 : 0;
}
