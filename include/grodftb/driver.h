/*
 * SDD:specs.md:§5.3 — Public C API for the DFTB+ driver library.
 *
 * Lifecycle: grodftb_init() → [set_geometry → compute → get_*]* → grodftb_finalize()
 */

#ifndef GRODFTB_DRIVER_H
#define GRODFTB_DRIVER_H

#include "grodftb/error.h"

#ifdef __cplusplus
extern "C" {
#endif

/* SDD:specs.md:§5.3 — grodftb_handle_t typedef is in error.h (needed by error API) */

/* SDD:specs.md:§17 — Version info struct */
typedef struct {
    int major;
    int minor;
    int patch;
    const char *git_hash;
} grodftb_version_t;

/**
 * Get library version information.
 *
 * Returns a pointer to a static struct containing the version number
 * and git commit hash. The pointer is valid for the lifetime of the
 * program and must not be freed by the caller.
 *
 * SDD:specs.md:§5.3
 *
 * @return Pointer to static version struct (never NULL)
 */
const grodftb_version_t* grodftb_version(void);

/**
 * Initialize a DFTB+ instance from an HSD template file.
 *
 * Executes the 4-step DFTB+ init sequence:
 *   dftbp_init → dftbp_get_input_from_file → dftbp_process_input → dftbp_input_final
 *
 * The HSD template defines the DFTB method, SK parameters, SCC settings,
 * and species list. Geometry and embedding are set separately via
 * grodftb_set_geometry() and grodftb_set_embedding_*().
 *
 * @param template_path  Path to .hsd template file
 * @param handle_out     Receives the new handle on success
 * @return GRODFTB_SUCCESS on success, error code otherwise
 */
int grodftb_init(const char *template_path, grodftb_handle_t *handle_out);

/**
 * Release all resources associated with the handle.
 * Safe to call with NULL or pointer to NULL. Sets *handle to NULL.
 *
 * @param handle  Pointer to handle (set to NULL on return)
 */
void grodftb_finalize(grodftb_handle_t *handle);

/**
 * Set QM atom coordinates and species in the driver context.
 *
 * Copies coordinates and species into internal buffers, then propagates
 * coordinates to DFTB+ via dftbp_set_coords(). Must be called before
 * each compute(). Coordinates are in Bohr, flat layout [x0,y0,z0,...].
 *
 * ThDD:US-011:Eq1.1, Eq1.2, Eq1.4
 * SDD:specs.md:§5.3
 *
 * @param handle   Active driver handle from grodftb_init()
 * @param natoms   Number of QM atoms (must match init geometry)
 * @param species  Species indices, 0-based [natoms]
 * @param coords   Coordinates in Bohr, flat [3*natoms]
 * @return GRODFTB_SUCCESS on success, error code otherwise
 */
int grodftb_set_geometry(grodftb_handle_t handle, int natoms,
                         const int *species, const double *coords);

/**
 * Run the SCC-DFTB calculation to self-consistency.
 *
 * Triggers the full DFTB+ SCC solve at the current geometry, then
 * retrieves and caches energy, forces (negated gradients), and Mulliken
 * charges. All subsequent get_*() calls return data from this solve.
 *
 * Requires: grodftb_init() and grodftb_set_geometry() called prior.
 *
 * ThDD:US-012:Eq3.2 — compute sequence
 * ThDD:06_theory:Eq3.2 — DFTB2 energy functional
 * SDD:specs.md:§5.3
 *
 * @param handle  Active driver handle
 * @return GRODFTB_SUCCESS on success,
 *         GRODFTB_ERR_NULL_HANDLE if handle is NULL,
 *         GRODFTB_ERR_NOT_INITIALIZED if init was not called,
 *         GRODFTB_ERR_INVALID_ARGUMENT if geometry was not set,
 *         GRODFTB_ERR_SCC_NOT_CONVERGED if the SCC cycle fails,
 *         GRODFTB_ERR_DFTB_INTERNAL for other DFTB+ errors
 */
int grodftb_compute(grodftb_handle_t handle);

/**
 * Get total QM energy (Hartree) from the last grodftb_compute() call.
 *
 * ThDD:06_theory:Eq3.2 — Mermin free energy
 * SDD:specs.md:§5.3
 *
 * @param handle      Active driver handle
 * @param energy_out  Receives the Mermin free energy in Hartree
 * @return GRODFTB_SUCCESS on success,
 *         GRODFTB_ERR_NULL_HANDLE if handle is NULL,
 *         GRODFTB_ERR_NOT_INITIALIZED if init was not called,
 *         GRODFTB_ERR_NO_RESULTS if compute has not been called
 */
int grodftb_get_energy(grodftb_handle_t handle, double *energy_out);

/**
 * Get QM forces (Hartree/Bohr) from the last grodftb_compute() call.
 *
 * Forces are the negated DFTB+ gradients: F = -dE/dR.
 * Output is a flat array [x0,y0,z0,x1,y1,z1,...] of size 3*natoms.
 *
 * ThDD:06_theory:Eq4.1 — F_K = -∂E/∂R_K
 * SDD:specs.md:§5.3
 *
 * @param handle      Active driver handle
 * @param forces_out  Receives forces in Hartree/Bohr, flat [3*natoms]
 * @return GRODFTB_SUCCESS on success,
 *         GRODFTB_ERR_NULL_HANDLE if handle is NULL,
 *         GRODFTB_ERR_INVALID_ARGUMENT if forces_out is NULL,
 *         GRODFTB_ERR_NO_RESULTS if compute has not been called
 */
int grodftb_get_forces(grodftb_handle_t handle, double *forces_out);

/**
 * Retrieve Mulliken gross charges on QM atoms after compute.
 *
 * Charges follow the standard DFTB convention: negative values indicate
 * electron excess. Output units: elementary charge (e). No conversion needed.
 *
 * ThDD:06_theory:§2.3 — q_A = sum_{mu in A} sum_nu P_{mu,nu} S_{mu,nu}
 * SDD:specs.md:§5.3
 *
 * @param handle       Active driver handle
 * @param charges_out  Receives Mulliken charges [natoms], pre-allocated by caller
 * @return GRODFTB_SUCCESS on success,
 *         GRODFTB_ERR_NULL_HANDLE if handle is NULL,
 *         GRODFTB_ERR_INVALID_ARGUMENT if charges_out is NULL,
 *         GRODFTB_ERR_NO_RESULTS if compute has not been called
 */
int grodftb_get_mulliken_charges(grodftb_handle_t handle, double *charges_out);

/**
 * Check whether the last grodftb_compute() converged.
 *
 * ThDD:US-012:Eq2.6 — SCC convergence criterion
 * SDD:specs.md:§5.3
 *
 * @param handle         Active driver handle
 * @param converged_out  Receives 1 if converged, 0 otherwise
 * @return GRODFTB_SUCCESS, or GRODFTB_ERR_NULL_HANDLE / GRODFTB_ERR_NOT_INITIALIZED
 */
int grodftb_is_converged(grodftb_handle_t handle, int *converged_out);

/**
 * Get the number of SCC iterations from the last grodftb_compute() call.
 *
 * Note: DFTB+ C API v0.4.0 does not expose iteration count.
 * Returns 0 (unknown) until a future DFTB+ API extension.
 *
 * SDD:specs.md:§5.3
 *
 * @param handle     Active driver handle
 * @param niter_out  Receives the iteration count (0 if unknown)
 * @return GRODFTB_SUCCESS, or GRODFTB_ERR_NULL_HANDLE / GRODFTB_ERR_NOT_INITIALIZED
 */
int grodftb_get_scc_iterations(grodftb_handle_t handle, int *niter_out);

/**
 * Set MM point charges for electrostatic embedding.
 *
 * Computes the Coulomb potential and its gradient at each QM atom site
 * from the given MM point charges, then passes to DFTB+ via
 * dftbp_set_external_potential(). All quantities in atomic units.
 *
 * Sign conventions (from DFTB+ extchargepot.f90):
 *   extpot[A]      = -sum_J Q_J / |R_A - R_J|
 *   extpotgrad[A,a] = +sum_J Q_J * (R_A - R_J)_a / |R_A - R_J|^3
 *
 * Call with ncharges=0 to clear embedding and revert to gas-phase.
 * Requires grodftb_set_geometry() to have been called first.
 *
 * ThDD:06_theory:Eq1.1 — Coulomb potential from MM point charges
 * SDD:specs.md:§5.3
 *
 * @param handle     Active driver handle
 * @param ncharges   Number of MM point charges (0 to clear)
 * @param charges    MM charges in elementary charge (e) [ncharges]
 * @param positions  MM positions in Bohr, flat [3*ncharges]
 * @return GRODFTB_SUCCESS on success,
 *         GRODFTB_ERR_NULL_HANDLE if handle is NULL,
 *         GRODFTB_ERR_NOT_INITIALIZED if init was not called,
 *         GRODFTB_ERR_INVALID_ARGUMENT if geometry not set, ncharges < 0,
 *             or charges/positions NULL when ncharges > 0
 */
int grodftb_set_embedding_charges(grodftb_handle_t handle,
                                  int ncharges,
                                  const double *charges,
                                  const double *positions);

/**
 * Set external electrostatic potential (and optional field) at QM atom sites.
 * SDD:specs.md:§5.3
 *
 * PME pathway: caller pre-computes potential and electric field from MM
 * environment; this function passes them to DFTB+ via
 * dftbp_set_external_potential().
 *
 * @param handle     Opaque context handle
 * @param natoms     Number of QM atoms (must match set_geometry)
 * @param potential  Electrostatic potential at each QM atom, Hartree [natoms]
 * @param field      Electric field at each QM atom, Hartree/Bohr [3*natoms],
 *                   or NULL if field is unavailable (forces will lack embedding
 *                   contribution)
 * @return GRODFTB_SUCCESS on success,
 *         GRODFTB_ERR_NULL_HANDLE if handle is NULL,
 *         GRODFTB_ERR_NOT_INITIALIZED if init was not called,
 *         GRODFTB_ERR_INVALID_ARGUMENT if geometry not set, potential NULL,
 *             or natoms does not match set_geometry
 */
int grodftb_set_embedding_potential(grodftb_handle_t handle,
                                    int natoms,
                                    const double *potential,
                                    const double *field);

/**
 * Set embedding parameters for cutoff-based switching.
 *
 * Configures the cutoff radius and optional switching function width.
 * When switch_width > 0, a quintic switching function smoothly reduces
 * the embedding interaction to zero near the cutoff boundary, eliminating
 * force discontinuities that cause energy drift in NVE simulations.
 *
 * The switching function S(u) = 1 - 10u³ + 15u⁴ - 6u⁵ is applied where
 * u = (r - r_on) / switch_width, r_on = cutoff - switch_width.
 *
 * This must be called BEFORE grodftb_set_embedding_charges() for the
 * switching to take effect on the embedding potential and forces.
 *
 * ThDD:T-US-039-2.1 — Quintic switching function
 * ThDD:T-US-039-U.1 — Units: all parameters in Bohr
 * SDD:specs.md:§8.1
 *
 * @param handle       Active driver handle
 * @param cutoff_bohr  Cutoff radius in Bohr (must be > 0)
 * @param switch_width_bohr  Switching region width in Bohr.
 *                     Set to 0 for hard cutoff (no switching).
 *                     Must satisfy: 0 <= switch_width < cutoff
 * @return GRODFTB_SUCCESS on success,
 *         GRODFTB_ERR_NULL_HANDLE if handle is NULL,
 *         GRODFTB_ERR_NOT_INITIALIZED if init was not called,
 *         GRODFTB_ERR_INVALID_ARGUMENT if cutoff <= 0 or switch_width >= cutoff
 */
int grodftb_set_embedding_params(grodftb_handle_t handle,
                                 double cutoff_bohr,
                                 double switch_width_bohr);

/**
 * Set Gaussian damping width for electrostatic embedding.
 *
 * When sigma > 0, the bare Coulomb potential 1/r is replaced by the
 * Gaussian-damped erf(r/sigma)/r to eliminate short-range divergence
 * that causes SCC convergence failures when MM atoms approach QM atoms.
 *
 * sigma = 0 disables damping (bare Coulomb, default behavior).
 *
 * This must be called BEFORE grodftb_set_embedding_charges() for the
 * damping to take effect.
 *
 * ThDD:T-US-043b-6.3 — sigma=0 default guard
 * SDD:specs.md:§5.3
 *
 * @param handle      Active driver handle
 * @param sigma_bohr  Damping width in Bohr (must be >= 0)
 * @return GRODFTB_SUCCESS on success,
 *         GRODFTB_ERR_NULL_HANDLE if handle is NULL,
 *         GRODFTB_ERR_NOT_INITIALIZED if init was not called,
 *         GRODFTB_ERR_INVALID_ARGUMENT if sigma_bohr < 0
 */
int grodftb_set_embedding_damping(grodftb_handle_t handle,
                                   double sigma_bohr);

/**
 * Get Coulomb back-reaction forces on MM point charges (Hartree/Bohr).
 *
 * Computes F_{J,α} = +Q_J Σ_A Δq_A (R_J - R_A)_α / |R_J - R_A|³
 * from converged Mulliken charges and stored MM charge data.
 * Output is a flat array [x0,y0,z0,...] of size 3*ncharges.
 *
 * If switching is enabled via grodftb_set_embedding_params(), forces are
 * computed using the switched formula (ThDD:T-US-039-6.5):
 *   F_J = Q_J Σ_A q_A [S/r³ - S'/(w·r²)] (R_J - R_A)
 * which ensures smooth forces at the cutoff boundary.
 *
 * Only valid after grodftb_compute() when embedding charges were set
 * via grodftb_set_embedding_charges().
 *
 * ThDD:T-US-021-3.6 — Coulomb force on MM charges (unswitched)
 * ThDD:T-US-039-6.5 — Switched force formula
 * ThDD:T-US-021-3.3 — Hellmann-Feynman: no Pulay terms at SCC convergence
 * SDD:specs.md:§5.3
 *
 * @param handle      Active driver handle
 * @param forces_out  Receives forces in Hartree/Bohr, flat [3*ncharges]
 * @return GRODFTB_SUCCESS on success,
 *         GRODFTB_ERR_NULL_HANDLE if handle is NULL,
 *         GRODFTB_ERR_NOT_INITIALIZED if init was not called,
 *         GRODFTB_ERR_INVALID_ARGUMENT if forces_out is NULL or no embedding set,
 *         GRODFTB_ERR_NO_RESULTS if compute has not been called
 */
int grodftb_get_embedding_forces(grodftb_handle_t handle, double *forces_out);

#ifdef __cplusplus
}
#endif

#endif /* GRODFTB_DRIVER_H */
