/*
 * ThDD:US-009:Eq1.1–1.5 — DFTB+ driver init/finalize implementation.
 * ThDD:US-019 — Electrostatic embedding via dftbp_set_external_potential().
 *
 * Wraps the 4-step DFTB+ C API initialization sequence:
 *   1. dftbp_init()                 — create instance (Eq1.1)
 *   2. dftbp_get_input_from_file()  — parse HSD (Eq1.2)
 *   3. dftbp_process_input()        — set up calculator (Eq1.3)
 *   4. dftbp_input_final()          — release input tree (Eq1.4)
 *
 * ThDD:US-009:Eq3.2 — grodftb_finalize() guarantees all resources freed.
 *
 * SDD:specs.md:§5.3 — Function signatures match specification.
 * SDD:specs.md:§5.4 — Lifecycle contract enforced.
 * SDD:specs.md:§5.5 — Memory management: finalize frees all.
 */

#include "grodftb/driver.h"
#include "grodftb/damping.h"
#include "grodftb/switching.h"

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef GRODFTB_HAS_DFTBPLUS
#include <dftbplus.h>
#endif

/*
 * SDD:specs.md:§5.3 (context) — Internal state for a DFTB+ instance.
 *
 * ThDD:US-009:Eq2.1 — State machine: UNINITIALIZED → INITIALIZED → ...
 */
struct grodftb_context {
#ifdef GRODFTB_HAS_DFTBPLUS
    DftbPlus dftb_instance;
#endif
    int natoms;
    int ncharges;
    int *species;
    double *coords;
    double lattice[9];
    int periodic;
    double *ext_charges;
    double *ext_positions;
    double energy;
    double *forces;
    double *emb_forces;
    double *ext_potential;
    double *ext_potential_grad;
    double *mulliken;
    int scc_converged;
    int scc_iterations;
    int initialized;
    char last_error_msg[256];
    /*
     * ThDD:T-US-039-N.13 — Switching function parameters for cutoff embedding
     * Units: Bohr. When switch_width_bohr = 0, hard cutoff is used.
     */
    double cutoff_bohr;
    double switch_width_bohr;
    double r_on_bohr;           /* = cutoff_bohr - switch_width_bohr */
    /*
     * ThDD:T-US-043b-6.3 — Gaussian damping parameters for embedding
     * Units: Bohr. When damp_sigma = 0, bare Coulomb is used.
     * Zero-initialized by calloc (use_damping = 0, damp_sigma = 0.0).
     */
    double damp_sigma;
    int use_damping;
};

/* Check if a file exists and is readable. */
static int file_exists(const char *path)
{
    FILE *f = fopen(path, "r");
    if (f) {
        fclose(f);
        return 1;
    }
    return 0;
}

/*
 * SDD:specs.md:§18.1 — Set per-handle error message with printf-style formatting.
 * Thread-local: stored in ctx->last_error_msg, valid until next grodftb_* call.
 */
static void set_error(struct grodftb_context *ctx, const char *fmt, ...)
{
    if (!ctx) return;
    va_list args;
    va_start(args, fmt);
    vsnprintf(ctx->last_error_msg, sizeof(ctx->last_error_msg), fmt, args);
    va_end(args);
}

/*
 * SDD:specs.md:§5.3 — grodftb_init()
 * ThDD:US-009:Eq1.1–1.5 — 4-step DFTB+ initialization sequence.
 */
int grodftb_init(const char *template_path, grodftb_handle_t *handle_out)
{
    struct grodftb_context *ctx;

    /* Validate arguments (ThDD:US-009:Eq4.1) */
    if (!template_path || !handle_out) {
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }

    *handle_out = NULL;

    /* Check file exists before calling DFTB+ (which may abort on bad input) */
    if (!file_exists(template_path)) {
        return GRODFTB_ERR_FILE_NOT_FOUND;
    }

    /* Allocate context */
    ctx = calloc(1, sizeof(*ctx));
    if (!ctx) {
        return GRODFTB_ERR_ALLOC_FAILED;
    }

#ifdef GRODFTB_HAS_DFTBPLUS
    /* ThDD:US-009:Eq1.1 — Create DFTB+ instance.
     * Second argument is the output filename, not the HSD path.
     * Use "dftb_debug.out" for debugging, "/dev/null" to suppress. */
    dftbp_init(&ctx->dftb_instance, "/dev/null");

    /* ThDD:US-009:Eq1.2 — Parse HSD file into input tree. */
    {
        DftbPlusInput input;
        dftbp_get_input_from_file(&ctx->dftb_instance, template_path, &input);

        /* ThDD:US-009:Eq1.3 — Set up calculator (allocates matrices, reads SK files). */
        dftbp_process_input(&ctx->dftb_instance, &input);

        /* ThDD:US-009:Eq1.4 — Release input tree. */
        dftbp_input_final(&input);
    }

    ctx->natoms = dftbp_get_nr_atoms(&ctx->dftb_instance);
#endif

    ctx->initialized = 1;
    *handle_out = ctx;
    return GRODFTB_SUCCESS;
}

/*
 * SDD:specs.md:§5.3 — grodftb_finalize()
 * ThDD:US-009:Eq3.2 — Resource cleanup guarantee.
 * Safe to call with NULL or pointer-to-NULL.
 */
void grodftb_finalize(grodftb_handle_t *handle)
{
    struct grodftb_context *ctx;

    if (!handle || !*handle) {
        return;
    }

    ctx = *handle;

#ifdef GRODFTB_HAS_DFTBPLUS
    if (ctx->initialized) {
        dftbp_final(&ctx->dftb_instance);
    }
#endif

    free(ctx->species);
    free(ctx->coords);
    free(ctx->ext_charges);
    free(ctx->ext_positions);
    free(ctx->forces);
    free(ctx->emb_forces);
    free(ctx->ext_potential);
    free(ctx->ext_potential_grad);
    free(ctx->mulliken);

    free(ctx);
    *handle = NULL;
}

/*
 * SDD:specs.md:§5.3 — grodftb_set_geometry()
 * ThDD:US-011:Eq1.1 — Flat coordinate layout [x0,y0,z0,x1,y1,z1,...]
 * ThDD:US-011:Eq1.2 — Fortran (3,natom) = C [natom][3] layout equivalence
 * ThDD:US-011:Eq1.4 — Operation sequence: validate → copy → propagate → invalidate
 * SDD:specs.md:§5.5 — Allocate on first use, steady-state after
 */
int grodftb_set_geometry(grodftb_handle_t handle, int natoms,
                         const int *species, const double *coords)
{
    struct grodftb_context *ctx = handle;

    /* ThDD:US-011:Eq3.1 — Error condition checks */
    if (!ctx) {
        return GRODFTB_ERR_NULL_HANDLE;
    }
    ctx->last_error_msg[0] = '\0';
    if (!ctx->initialized) {
        set_error(ctx, "handle not initialized — call grodftb_init() first");
        return GRODFTB_ERR_NOT_INITIALIZED;
    }
    if (natoms <= 0) {
        set_error(ctx, "natoms=%d is invalid (must be > 0)", natoms);
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    if (!species) {
        set_error(ctx, "species pointer is NULL");
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    if (!coords) {
        set_error(ctx, "coords pointer is NULL");
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    if (natoms != ctx->natoms) {
        set_error(ctx, "natoms mismatch: expected %d, got %d", ctx->natoms, natoms);
        return GRODFTB_ERR_SIZE_MISMATCH;
    }

    /* SDD:specs.md:§5.5 — Allocate buffers on first call */
    if (!ctx->coords) {
        ctx->coords = malloc(3 * natoms * sizeof(double));
        if (!ctx->coords) {
            return GRODFTB_ERR_ALLOC_FAILED;
        }
    }
    if (!ctx->species) {
        ctx->species = malloc(natoms * sizeof(int));
        if (!ctx->species) {
            return GRODFTB_ERR_ALLOC_FAILED;
        }
    }
    if (!ctx->forces) {
        ctx->forces = malloc(3 * natoms * sizeof(double));
        if (!ctx->forces) {
            return GRODFTB_ERR_ALLOC_FAILED;
        }
    }
    if (!ctx->mulliken) {
        ctx->mulliken = malloc(natoms * sizeof(double));
        if (!ctx->mulliken) {
            return GRODFTB_ERR_ALLOC_FAILED;
        }
    }

    /* ThDD:US-011:Eq1.1 — Copy coordinates (bit-exact via memcpy) */
    memcpy(ctx->coords, coords, 3 * natoms * sizeof(double));

    /* ThDD:US-011:Eq1.4 — Copy species */
    memcpy(ctx->species, species, natoms * sizeof(int));

    /* ThDD:US-011:Eq1.4 — Propagate coordinates to DFTB+ */
#ifdef GRODFTB_HAS_DFTBPLUS
    dftbp_set_coords(&ctx->dftb_instance, ctx->coords);
#endif

    /* Invalidate cached results — force fresh SCC solve on next compute() */
    ctx->scc_converged = 0;

    return GRODFTB_SUCCESS;
}

/*
 * SDD:specs.md:§5.3 — grodftb_compute()
 * ThDD:US-012:Eq3.2 — Compute sequence: energy → gradients → negate → charges
 * ThDD:06_theory:Eq3.2 — DFTB2 total energy (Mermin free energy)
 * ThDD:06_theory:Eq4.1 — F_A = -dE/dR_A (negate gradients to get forces)
 *
 * No heap allocation: all buffers allocated by grodftb_set_geometry().
 * dftbp_get_energy() triggers the SCC cycle internally (T-US-012-3.1).
 */
int grodftb_compute(grodftb_handle_t handle)
{
    struct grodftb_context *ctx = handle;
    int i;

    if (!ctx) {
        return GRODFTB_ERR_NULL_HANDLE;
    }
    ctx->last_error_msg[0] = '\0';
    if (!ctx->initialized) {
        set_error(ctx, "handle not initialized");
        return GRODFTB_ERR_NOT_INITIALIZED;
    }
    if (!ctx->coords) {
        set_error(ctx, "geometry not set — call grodftb_set_geometry() first");
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }

    /* Reset results from previous compute */
    ctx->scc_converged = 0;
    ctx->scc_iterations = 0;

#ifdef GRODFTB_HAS_DFTBPLUS
    /* ThDD:06_theory:Eq3.2 — Trigger SCC solve; retrieve Mermin free energy */
    dftbp_get_energy(&ctx->dftb_instance, &ctx->energy);

    /* ThDD:06_theory:Eq4.1 — Retrieve gradients (+dE/dR) */
    dftbp_get_gradients(&ctx->dftb_instance, ctx->forces);

    /* ThDD:06_theory:Eq4.1 — Negate gradients to obtain forces: F = -dE/dR */
    for (i = 0; i < 3 * ctx->natoms; i++) {
        ctx->forces[i] = -ctx->forces[i];
    }

    /* ThDD:US-011:Eq2.1 — Retrieve Mulliken gross charges */
    dftbp_get_gross_charges(&ctx->dftb_instance, ctx->mulliken);

    /* ThDD:T-US-021-3.6 — Embedding forces on MM charges (unswitched)
     * ThDD:T-US-039-6.5 — Switched embedding forces on MM charges
     *
     * Unswitched: F_J = Q_J Σ_A Δq_A (R_J - R_A) / r³
     * Switched:   F_J = Q_J Σ_A q_A [S/r³ - S'/(w·r²)] (R_J - R_A)
     *
     * ThDD:T-US-021-3.3 — Hellmann-Feynman theorem: at SCC convergence,
     * Pulay terms vanish, so analytic Coulomb formula is exact.
     */
    if (ctx->ncharges > 0 && ctx->emb_forces) {
        int j, a;
        const int use_switching = (ctx->switch_width_bohr > 0.0);
        const double r_on = ctx->r_on_bohr;
        const double w = ctx->switch_width_bohr;
        const double cutoff = ctx->cutoff_bohr;

        for (j = 0; j < ctx->ncharges; j++) {
            double fj[3] = {0.0, 0.0, 0.0};
            for (a = 0; a < ctx->natoms; a++) {
                /* ThDD:T-US-021-3.5 — Coulomb kernel derivative */
                double dx = ctx->ext_positions[3*j + 0] - ctx->coords[3*a + 0];
                double dy = ctx->ext_positions[3*j + 1] - ctx->coords[3*a + 1];
                double dz = ctx->ext_positions[3*j + 2] - ctx->coords[3*a + 2];
                double r2 = dx*dx + dy*dy + dz*dz;
                double r  = sqrt(r2);

                if (use_switching) {
                    /*
                     * ThDD:T-US-039-6.5 — Switched force formula
                     * ThDD:T-US-043b-4.6 — Damped switched force formula
                     */
                    if (r >= cutoff) {
                        continue;
                    }

                    double S, dSdr;
                    grodftb_switch_func_and_deriv(r, r_on, w, &S, &dSdr);

                    double kernel;
                    if (ctx->use_damping) {
                        /* ThDD:T-US-043b-4.6 — Damped+switched force kernel */
                        double f_val, K_damp;
                        grodftb_damped_coulomb_and_kernel(r, ctx->damp_sigma,
                                                          &f_val, &K_damp);
                        if (r < 1e-10) {
                            kernel = S * K_damp;
                        } else {
                            kernel = S * K_damp - dSdr * f_val / r;
                        }
                    } else {
                        /* ThDD:T-US-039-6.5 — Bare Coulomb switched kernel */
                        double inv_r = 1.0 / r;
                        double inv_r2 = inv_r * inv_r;
                        double inv_r3 = inv_r2 * inv_r;
                        kernel = S * inv_r3 - dSdr * inv_r2;
                    }

                    double q_kernel = ctx->mulliken[a] * kernel;
                    fj[0] += q_kernel * dx;
                    fj[1] += q_kernel * dy;
                    fj[2] += q_kernel * dz;
                } else {
                    double kernel;
                    if (ctx->use_damping) {
                        /* ThDD:T-US-043b-4.6 — Unswitched damped force */
                        double f_val, K_damp;
                        grodftb_damped_coulomb_and_kernel(r, ctx->damp_sigma,
                                                          &f_val, &K_damp);
                        kernel = K_damp;
                    } else {
                        /* ThDD:T-US-021-3.6 — Unswitched Coulomb force */
                        kernel = 1.0 / (r * r2);
                    }
                    double dq_kernel = ctx->mulliken[a] * kernel;
                    fj[0] += dq_kernel * dx;
                    fj[1] += dq_kernel * dy;
                    fj[2] += dq_kernel * dz;
                }
            }
            /* ThDD:T-US-021-3.6 — multiply by +Q_J
             * E_emb = -Σ_{A,J} Δq_A Q_J / r_AJ (from Eq T-US-021-2.3)
             * dE/dR_J = -Σ_A Δq_A Q_J (R_J - R_A) / r^3
             * F_J = -dE/dR_J = +Q_J Σ_A Δq_A (R_J - R_A) / r^3
             * FD-validated: relative error < 10^-9 on B4. */
            ctx->emb_forces[3*j + 0] = ctx->ext_charges[j] * fj[0];
            ctx->emb_forces[3*j + 1] = ctx->ext_charges[j] * fj[1];
            ctx->emb_forces[3*j + 2] = ctx->ext_charges[j] * fj[2];
        }
    }

    ctx->scc_converged = 1;
#else
    /* Stub mode: zero all results */
    ctx->energy = 0.0;
    memset(ctx->forces, 0, 3 * ctx->natoms * sizeof(double));
    memset(ctx->mulliken, 0, ctx->natoms * sizeof(double));
    if (ctx->ncharges > 0 && ctx->emb_forces) {
        memset(ctx->emb_forces, 0, 3 * ctx->ncharges * sizeof(double));
    }
    ctx->scc_converged = 1;
#endif

    return GRODFTB_SUCCESS;
}

/*
 * SDD:specs.md:§5.3 — grodftb_get_energy()
 * ThDD:06_theory:Eq3.2 — Returns Mermin free energy in Hartree.
 */
int grodftb_get_energy(grodftb_handle_t handle, double *energy_out)
{
    struct grodftb_context *ctx = handle;

    if (!ctx) {
        return GRODFTB_ERR_NULL_HANDLE;
    }
    ctx->last_error_msg[0] = '\0';
    if (!ctx->initialized) {
        set_error(ctx, "handle not initialized");
        return GRODFTB_ERR_NOT_INITIALIZED;
    }
    if (!energy_out) {
        set_error(ctx, "output pointer is NULL");
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    if (!ctx->scc_converged) {
        set_error(ctx, "no results — call grodftb_compute() first");
        return GRODFTB_ERR_NO_RESULTS;
    }

    *energy_out = ctx->energy;
    return GRODFTB_SUCCESS;
}

/*
 * SDD:specs.md:§5.3 — grodftb_get_forces()
 * ThDD:06_theory:Eq4.1 — F_K = -dE/dR_K (negated gradients, cached by compute)
 * ThDD:T-US-014-5.1 — Pass-through: memcpy from ctx->forces
 *
 * Forces are in Hartree/Bohr, flat [3*natoms]. No unit conversion.
 */
int grodftb_get_forces(grodftb_handle_t handle, double *forces_out)
{
    struct grodftb_context *ctx = handle;

    if (!ctx) {
        return GRODFTB_ERR_NULL_HANDLE;
    }
    ctx->last_error_msg[0] = '\0';
    if (!ctx->initialized) {
        set_error(ctx, "handle not initialized");
        return GRODFTB_ERR_NOT_INITIALIZED;
    }
    if (!forces_out) {
        set_error(ctx, "output pointer is NULL");
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    if (!ctx->scc_converged) {
        set_error(ctx, "no results — call grodftb_compute() first");
        return GRODFTB_ERR_NO_RESULTS;
    }

    memcpy(forces_out, ctx->forces, 3 * ctx->natoms * sizeof(double));
    return GRODFTB_SUCCESS;
}

/*
 * SDD:specs.md:§5.3 — grodftb_get_mulliken_charges()
 * ThDD:06_theory:§2.3 — Mulliken gross charges q_A = sum_{mu in A} (PS)_{mu,mu}
 * ThDD:T-US-016-2.1 — Pass-through: memcpy from ctx->mulliken
 *
 * Charges are in elementary charge (e), flat [natoms]. No unit conversion.
 * Sign convention: negative = electron excess (standard DFTB).
 */
int grodftb_get_mulliken_charges(grodftb_handle_t handle, double *charges_out)
{
    struct grodftb_context *ctx = handle;

    if (!ctx) {
        return GRODFTB_ERR_NULL_HANDLE;
    }
    ctx->last_error_msg[0] = '\0';
    if (!ctx->initialized) {
        set_error(ctx, "handle not initialized");
        return GRODFTB_ERR_NOT_INITIALIZED;
    }
    if (!charges_out) {
        set_error(ctx, "output pointer is NULL");
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    if (!ctx->scc_converged) {
        set_error(ctx, "no results — call grodftb_compute() first");
        return GRODFTB_ERR_NO_RESULTS;
    }

    memcpy(charges_out, ctx->mulliken, ctx->natoms * sizeof(double));
    return GRODFTB_SUCCESS;
}

/*
 * SDD:specs.md:§8.1 — grodftb_set_embedding_params()
 * ThDD:T-US-039-U.1 — Units: all parameters in Bohr
 * ThDD:T-US-039-U.2 — r_on = cutoff - switch_width
 * ThDD:T-US-039-N.13 — Storage in context
 *
 * Configures cutoff-based embedding with optional quintic switching.
 */
int grodftb_set_embedding_params(grodftb_handle_t handle,
                                 double cutoff_bohr,
                                 double switch_width_bohr)
{
    struct grodftb_context *ctx = handle;

    if (!ctx)
        return GRODFTB_ERR_NULL_HANDLE;
    ctx->last_error_msg[0] = '\0';
    if (!ctx->initialized) {
        set_error(ctx, "handle not initialized");
        return GRODFTB_ERR_NOT_INITIALIZED;
    }

    /* ThDD:T-US-039-U.1 — Validate parameters */
    if (cutoff_bohr <= 0.0) {
        set_error(ctx, "cutoff must be positive, got %g Bohr", cutoff_bohr);
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    if (switch_width_bohr < 0.0) {
        set_error(ctx, "switch_width must be non-negative, got %g Bohr", switch_width_bohr);
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    if (switch_width_bohr >= cutoff_bohr) {
        set_error(ctx, "switch_width (%g) must be less than cutoff (%g)",
                  switch_width_bohr, cutoff_bohr);
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }

    /* ThDD:T-US-039-U.2 — Compute r_on */
    ctx->cutoff_bohr = cutoff_bohr;
    ctx->switch_width_bohr = switch_width_bohr;
    ctx->r_on_bohr = cutoff_bohr - switch_width_bohr;

    return GRODFTB_SUCCESS;
}

/*
 * SDD:specs.md:§5.3 — grodftb_set_embedding_damping()
 * ThDD:T-US-043b-6.3 — sigma=0 disables damping (bare Coulomb default)
 *
 * Configures Gaussian damping width for electrostatic embedding.
 * When sigma > 0, erf(r/sigma)/r replaces 1/r in embedding potential.
 */
int grodftb_set_embedding_damping(grodftb_handle_t handle,
                                   double sigma_bohr)
{
    struct grodftb_context *ctx = handle;

    if (!ctx)
        return GRODFTB_ERR_NULL_HANDLE;
    ctx->last_error_msg[0] = '\0';
    if (!ctx->initialized) {
        set_error(ctx, "handle not initialized");
        return GRODFTB_ERR_NOT_INITIALIZED;
    }

    /* ThDD:T-US-043b-6.3 — Validate sigma */
    if (sigma_bohr < 0.0) {
        set_error(ctx, "sigma must be non-negative, got %g Bohr", sigma_bohr);
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }

    ctx->damp_sigma = sigma_bohr;
    ctx->use_damping = (sigma_bohr > 0.0);

    return GRODFTB_SUCCESS;
}

/*
 * SDD:specs.md:§5.3 — grodftb_set_embedding_charges()
 * ThDD:06_theory:Eq1.1 — Coulomb potential from MM point charges
 * ThDD:T-US-019-2.1 — extpot[A] = -sum_J Q_J / |R_A - R_J|
 * ThDD:T-US-019-3.1 — extpotgrad[A,a] = +sum_J Q_J * (R_A - R_J)_a / |R_A - R_J|^3
 *
 * Sign conventions match DFTB+ extchargepot.f90:61-62.
 */
int grodftb_set_embedding_charges(grodftb_handle_t handle,
                                  int ncharges,
                                  const double *charges,
                                  const double *positions)
{
    struct grodftb_context *ctx = handle;
    int i, k;

    if (!ctx) {
        return GRODFTB_ERR_NULL_HANDLE;
    }
    ctx->last_error_msg[0] = '\0';
    if (!ctx->initialized) {
        set_error(ctx, "handle not initialized — call grodftb_init() first");
        return GRODFTB_ERR_NOT_INITIALIZED;
    }
    if (!ctx->coords) {
        set_error(ctx, "geometry not set — call grodftb_set_geometry() first");
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    if (ncharges < 0) {
        set_error(ctx, "ncharges=%d is invalid (must be >= 0)", ncharges);
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    if (ncharges > 0 && !charges) {
        set_error(ctx, "charges pointer is NULL");
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    if (ncharges > 0 && !positions) {
        set_error(ctx, "positions pointer is NULL");
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }

    /* Ensure ext_potential and ext_potential_grad are allocated */
    if (!ctx->ext_potential) {
        ctx->ext_potential = calloc(ctx->natoms, sizeof(double));
        if (!ctx->ext_potential) {
            return GRODFTB_ERR_ALLOC_FAILED;
        }
    }
    if (!ctx->ext_potential_grad) {
        ctx->ext_potential_grad = calloc(3 * ctx->natoms, sizeof(double));
        if (!ctx->ext_potential_grad) {
            return GRODFTB_ERR_ALLOC_FAILED;
        }
    }

    /* Handle ncharges == 0: clear embedding */
    if (ncharges == 0) {
        memset(ctx->ext_potential, 0, ctx->natoms * sizeof(double));
        memset(ctx->ext_potential_grad, 0, 3 * ctx->natoms * sizeof(double));
        free(ctx->emb_forces);
        ctx->emb_forces = NULL;
        ctx->ncharges = 0;
#ifdef GRODFTB_HAS_DFTBPLUS
        dftbp_set_external_potential(&ctx->dftb_instance,
                                     ctx->ext_potential,
                                     ctx->ext_potential_grad);
#endif
        ctx->scc_converged = 0;
        return GRODFTB_SUCCESS;
    }

    /* Allocate/realloc charge storage if needed */
    if (ncharges != ctx->ncharges) {
        free(ctx->ext_charges);
        ctx->ext_charges = malloc(ncharges * sizeof(double));
        if (!ctx->ext_charges) {
            return GRODFTB_ERR_ALLOC_FAILED;
        }
        free(ctx->ext_positions);
        ctx->ext_positions = malloc(3 * ncharges * sizeof(double));
        if (!ctx->ext_positions) {
            return GRODFTB_ERR_ALLOC_FAILED;
        }
        free(ctx->emb_forces);
        ctx->emb_forces = malloc(3 * ncharges * sizeof(double));
        if (!ctx->emb_forces) {
            return GRODFTB_ERR_ALLOC_FAILED;
        }
    }

    /* Copy charges and positions */
    memcpy(ctx->ext_charges, charges, ncharges * sizeof(double));
    memcpy(ctx->ext_positions, positions, 3 * ncharges * sizeof(double));
    ctx->ncharges = ncharges;

    /* Zero accumulators */
    memset(ctx->ext_potential, 0, ctx->natoms * sizeof(double));
    memset(ctx->ext_potential_grad, 0, 3 * ctx->natoms * sizeof(double));

    /* ThDD:06_theory:Eq1.1 — Coulomb potential from MM point charges (unswitched)
     * ThDD:T-US-039-5.1 — Switched potential: V_A = -Σ_J Q_J S(r_AJ) / r_AJ
     * ThDD:T-US-039-5.2 — Switched gradient: ∂V_A/∂R_A = +Σ_J Q_J [S/r³ - S'/r²](R_A - R_J)
     */
    {
        const int use_switching = (ctx->switch_width_bohr > 0.0);
        const double r_on = ctx->r_on_bohr;
        const double w = ctx->switch_width_bohr;
        const double cutoff = ctx->cutoff_bohr;

        for (i = 0; i < ctx->natoms; i++) {
            for (k = 0; k < ncharges; k++) {
                double dx = ctx->coords[3*i + 0] - positions[3*k + 0];
                double dy = ctx->coords[3*i + 1] - positions[3*k + 1];
                double dz = ctx->coords[3*i + 2] - positions[3*k + 2];
                double r2 = dx*dx + dy*dy + dz*dz;
                double r  = sqrt(r2);

                /* ThDD:T-US-043b-9.1 — Overlap guard: only for bare Coulomb.
                 * When damping is active, erf(r/sigma)/r is finite at r=0
                 * (Taylor expansion handles it). */
                if (!ctx->use_damping && r < 1e-10) {
                    set_error(ctx, "QM atom %d and MM charge %d overlap (r=%.3e Bohr)",
                              i, k, r);
                    return GRODFTB_ERR_INVALID_ARGUMENT;
                }

                if (use_switching) {
                    /*
                     * ThDD:T-US-039-N.3 — Three-region handling for potential
                     * - r <= r_on: Full contribution, S=1
                     * - r_on < r < cutoff: Switched contribution
                     * - r >= cutoff: No contribution, S=0
                     */
                    if (r >= cutoff) {
                        /* Outside cutoff: no contribution */
                        continue;
                    }

                    double S, dSdr;
                    grodftb_switch_func_and_deriv(r, r_on, w, &S, &dSdr);

                    if (ctx->use_damping) {
                        /*
                         * ThDD:T-US-043b-3.5, ThDD:T-US-043b-3.6
                         * Damped+switched potential and gradient.
                         * kernel = -(d/dr[S*f])/r = S*K - dSdr*f/r
                         * where f=erf(r/sigma)/r, K=g/r, g=-df/dr.
                         */
                        double f_val, K_damp;
                        grodftb_damped_coulomb_and_kernel(r, ctx->damp_sigma,
                                                          &f_val, &K_damp);

                        /* ThDD:T-US-043b-3.5 — Damped+switched potential */
                        ctx->ext_potential[i] -= charges[k] * S * f_val;

                        /* ThDD:T-US-043b-3.6 — Damped+switched gradient kernel
                         * kernel = S * K_damp - dSdr * f_val / r
                         * where K_damp = g/r and dSdr is already dS/dr */
                        double kernel;
                        if (r < 1e-10) {
                            /* At r=0: dSdr=0 (since r < r_on), K_damp finite */
                            kernel = S * K_damp;
                        } else {
                            kernel = S * K_damp - dSdr * f_val / r;
                        }

                        ctx->ext_potential_grad[3*i + 0] += charges[k] * kernel * dx;
                        ctx->ext_potential_grad[3*i + 1] += charges[k] * kernel * dy;
                        ctx->ext_potential_grad[3*i + 2] += charges[k] * kernel * dz;
                    } else {
                        double inv_r = 1.0 / r;
                        double inv_r2 = inv_r * inv_r;
                        double inv_r3 = inv_r2 * inv_r;

                        /*
                         * ThDD:T-US-039-5.1 — Switched potential
                         * extpot[A] -= Q_J * S / r
                         */
                        ctx->ext_potential[i] -= charges[k] * S * inv_r;

                        /*
                         * ThDD:T-US-039-5.2 — Switched potential gradient
                         * kernel = S/r³ - dSdr/r²
                         * extpotgrad[A,a] += Q_J * kernel * (R_A - R_J)_a
                         */
                        double kernel = S * inv_r3 - dSdr * inv_r2;
                        ctx->ext_potential_grad[3*i + 0] += charges[k] * kernel * dx;
                        ctx->ext_potential_grad[3*i + 1] += charges[k] * kernel * dy;
                        ctx->ext_potential_grad[3*i + 2] += charges[k] * kernel * dz;
                    }
                } else {
                    if (ctx->use_damping) {
                        /* ThDD:T-US-043b-1.11 — Unswitched damped Coulomb */
                        double f_val, K_damp;
                        grodftb_damped_coulomb_and_kernel(r, ctx->damp_sigma,
                                                          &f_val, &K_damp);

                        ctx->ext_potential[i] -= charges[k] * f_val;

                        ctx->ext_potential_grad[3*i + 0] += charges[k] * K_damp * dx;
                        ctx->ext_potential_grad[3*i + 1] += charges[k] * K_damp * dy;
                        ctx->ext_potential_grad[3*i + 2] += charges[k] * K_damp * dz;
                    } else {
                        /* ThDD:T-US-019-2.1, T-US-019-3.1 — Unswitched Coulomb */
                        double inv_r = 1.0 / r;
                        double inv_r3 = inv_r / r2;

                        /* MINUS sign: extpot[A] -= Q_J / |R_A - R_J| */
                        ctx->ext_potential[i] -= charges[k] * inv_r;

                        /* PLUS sign: extpotgrad[A,a] += Q_J * (R_A - R_J)_a / |R_A - R_J|^3 */
                        ctx->ext_potential_grad[3*i + 0] += charges[k] * dx * inv_r3;
                        ctx->ext_potential_grad[3*i + 1] += charges[k] * dy * inv_r3;
                        ctx->ext_potential_grad[3*i + 2] += charges[k] * dz * inv_r3;
                    }
                }
            }
        }
    }

    /* Pass to DFTB+ */
#ifdef GRODFTB_HAS_DFTBPLUS
    dftbp_set_external_potential(&ctx->dftb_instance,
                                 ctx->ext_potential,
                                 ctx->ext_potential_grad);
#endif

    /* Invalidate cached results */
    ctx->scc_converged = 0;

    return GRODFTB_SUCCESS;
}

/*
 * SDD:specs.md:§5.3 — grodftb_set_embedding_potential()
 * ThDD:T-US-020-3.1 — extpot[A] = potential[A] (direct copy, no sign change)
 * ThDD:T-US-020-3.2 — extpotgrad[A,a] = -field[A,a] (sign negation: E = -∇φ)
 *
 * PME pathway: caller provides pre-computed potential and electric field at QM
 * atom sites. This function copies them into DFTB+'s extpot/extpotgrad buffers
 * and calls dftbp_set_external_potential(). Buffers are shared with US-019's
 * set_embedding_charges (allocated on first use by either function).
 *
 * Units: Hartree (potential), Hartree/Bohr (field) — all atomic units, no conversion.
 */
int grodftb_set_embedding_potential(grodftb_handle_t handle,
                                    int natoms,
                                    const double *potential,
                                    const double *field)
{
    struct grodftb_context *ctx = handle;
    int i;

    if (!ctx) {
        return GRODFTB_ERR_NULL_HANDLE;
    }
    ctx->last_error_msg[0] = '\0';
    if (!ctx->initialized) {
        set_error(ctx, "handle not initialized — call grodftb_init() first");
        return GRODFTB_ERR_NOT_INITIALIZED;
    }
    if (!ctx->coords) {
        set_error(ctx, "geometry not set — call grodftb_set_geometry() first");
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    if (!potential) {
        set_error(ctx, "potential pointer is NULL");
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    if (natoms != ctx->natoms) {
        set_error(ctx, "natoms mismatch: expected %d, got %d", ctx->natoms, natoms);
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }

    /* SDD:specs.md:§5.5 — Allocate on first use; reuse buffers */
    if (!ctx->ext_potential) {
        ctx->ext_potential = calloc(ctx->natoms, sizeof(double));
        if (!ctx->ext_potential) {
            return GRODFTB_ERR_ALLOC_FAILED;
        }
    }
    if (!ctx->ext_potential_grad) {
        ctx->ext_potential_grad = calloc(3 * ctx->natoms, sizeof(double));
        if (!ctx->ext_potential_grad) {
            return GRODFTB_ERR_ALLOC_FAILED;
        }
    }

    /* ThDD:T-US-020-3.1 — extpot[A] = potential[A] (direct copy) */
    memcpy(ctx->ext_potential, potential, ctx->natoms * sizeof(double));

    /* ThDD:T-US-020-3.2 — extpotgrad[i] = -field[i] (sign negation)
     * DFTB+ expects ∇φ (gradient of potential); caller provides E = -∇φ (electric field).
     * IEEE 754 sign-bit flip is exact — no rounding error introduced. */
    if (field) {
        for (i = 0; i < 3 * ctx->natoms; i++) {
            ctx->ext_potential_grad[i] = -field[i];
        }
    }

    /* Pass to DFTB+ */
#ifdef GRODFTB_HAS_DFTBPLUS
    dftbp_set_external_potential(&ctx->dftb_instance,
                                 ctx->ext_potential,
                                 field ? ctx->ext_potential_grad : NULL);
#endif

    /* Invalidate cached results */
    ctx->scc_converged = 0;

    return GRODFTB_SUCCESS;
}

/*
 * SDD:specs.md:§5.3 — grodftb_get_embedding_forces()
 * ThDD:T-US-021-3.6 — Forces on MM charges from QM-MM Coulomb interaction
 *
 * Returns forces (Hartree/Bohr) on MM point charges, flat [3*ncharges].
 * Forces computed eagerly in grodftb_compute() from converged Mulliken charges.
 */
int grodftb_get_embedding_forces(grodftb_handle_t handle, double *forces_out)
{
    struct grodftb_context *ctx = handle;

    if (!ctx)
        return GRODFTB_ERR_NULL_HANDLE;
    ctx->last_error_msg[0] = '\0';
    if (!ctx->initialized) {
        set_error(ctx, "handle not initialized");
        return GRODFTB_ERR_NOT_INITIALIZED;
    }
    if (!forces_out) {
        set_error(ctx, "output pointer is NULL");
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }
    if (!ctx->scc_converged) {
        set_error(ctx, "no results — call grodftb_compute() first");
        return GRODFTB_ERR_NO_RESULTS;
    }
    if (ctx->ncharges == 0) {
        set_error(ctx, "no embedding charges set — call grodftb_set_embedding_charges() first");
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }

    memcpy(forces_out, ctx->emb_forces, 3 * ctx->ncharges * sizeof(double));
    return GRODFTB_SUCCESS;
}

/*
 * SDD:specs.md:§5.3 — grodftb_is_converged()
 * ThDD:US-012:Eq2.6 — SCC convergence status.
 */
int grodftb_is_converged(grodftb_handle_t handle, int *converged_out)
{
    struct grodftb_context *ctx = handle;

    if (!ctx) {
        return GRODFTB_ERR_NULL_HANDLE;
    }
    ctx->last_error_msg[0] = '\0';
    if (!ctx->initialized) {
        set_error(ctx, "handle not initialized");
        return GRODFTB_ERR_NOT_INITIALIZED;
    }
    if (!converged_out) {
        set_error(ctx, "output pointer is NULL");
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }

    *converged_out = ctx->scc_converged;
    return GRODFTB_SUCCESS;
}

/*
 * SDD:specs.md:§5.3 — grodftb_get_scc_iterations()
 * DFTB+ C API v0.4.0 does not expose iteration count;
 * always returns 0 (unknown).
 */
int grodftb_get_scc_iterations(grodftb_handle_t handle, int *niter_out)
{
    struct grodftb_context *ctx = handle;

    if (!ctx) {
        return GRODFTB_ERR_NULL_HANDLE;
    }
    ctx->last_error_msg[0] = '\0';
    if (!ctx->initialized) {
        set_error(ctx, "handle not initialized");
        return GRODFTB_ERR_NOT_INITIALIZED;
    }
    if (!niter_out) {
        set_error(ctx, "output pointer is NULL");
        return GRODFTB_ERR_INVALID_ARGUMENT;
    }

    *niter_out = ctx->scc_iterations;
    return GRODFTB_SUCCESS;
}

/*
 * SDD:specs.md:§18.1 — grodftb_error_string()
 * Returns a static string name for each error code.
 */
const char *grodftb_error_string(int errcode)
{
    switch (errcode) {
    case GRODFTB_SUCCESS:                    return "GRODFTB_SUCCESS";
    case GRODFTB_ERR_NULL_HANDLE:            return "GRODFTB_ERR_NULL_HANDLE";
    case GRODFTB_ERR_NOT_INITIALIZED:        return "GRODFTB_ERR_NOT_INITIALIZED";
    case GRODFTB_ERR_ALREADY_INIT:           return "GRODFTB_ERR_ALREADY_INIT";
    case GRODFTB_ERR_INVALID_ARGUMENT:       return "GRODFTB_ERR_INVALID_ARGUMENT";
    case GRODFTB_ERR_SIZE_MISMATCH:          return "GRODFTB_ERR_SIZE_MISMATCH";
    case GRODFTB_ERR_FILE_NOT_FOUND:         return "GRODFTB_ERR_FILE_NOT_FOUND";
    case GRODFTB_ERR_HSD_PARSE:              return "GRODFTB_ERR_HSD_PARSE";
    case GRODFTB_ERR_DFTB_INIT:              return "GRODFTB_ERR_DFTB_INIT";
    case GRODFTB_ERR_SCC_NOT_CONVERGED:      return "GRODFTB_ERR_SCC_NOT_CONVERGED";
    case GRODFTB_ERR_NO_RESULTS:             return "GRODFTB_ERR_NO_RESULTS";
    case GRODFTB_ERR_EXCITED_NOT_CONFIGURED: return "GRODFTB_ERR_EXCITED_NOT_CONFIGURED";
    case GRODFTB_ERR_EXCITED_FAILED:         return "GRODFTB_ERR_EXCITED_FAILED";
    case GRODFTB_ERR_NAC_UNAVAILABLE:        return "GRODFTB_ERR_NAC_UNAVAILABLE";
    case GRODFTB_ERR_CP_SOLVE_FAILED:        return "GRODFTB_ERR_CP_SOLVE_FAILED";
    case GRODFTB_ERR_ALLOC_FAILED:           return "GRODFTB_ERR_ALLOC_FAILED";
    case GRODFTB_ERR_DFTB_INTERNAL:          return "GRODFTB_ERR_DFTB_INTERNAL";
    default:                                 return "GRODFTB_ERR_UNKNOWN";
    }
}

/*
 * SDD:specs.md:§18.1 — grodftb_last_error()
 * Returns a human-readable message for the last error on this handle.
 * Currently returns a generic message; per-handle error tracking will
 * be added when later stories introduce complex error paths.
 */
const char *grodftb_last_error(grodftb_handle_t handle)
{
    struct grodftb_context *ctx = handle;
    if (!ctx) {
        return "NULL handle — no error detail available";
    }
    if (ctx->last_error_msg[0] == '\0') {
        return "no error";
    }
    return ctx->last_error_msg;
}
