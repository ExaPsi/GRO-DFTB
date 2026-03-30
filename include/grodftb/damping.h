/*
 * SDD:specs.md:§8.1 -- Gaussian-damped Coulomb functions for embedding
 * ThDD:T-US-043b-1.11 -- Damped embedding potential: erf(r/sigma)/r
 * ThDD:T-US-043b-2.4  -- Radial derivative of damped potential
 * ThDD:T-US-043b-9.1  -- Taylor expansion for small r/sigma
 *
 * These helper functions compute the Gaussian-damped Coulomb potential
 * erf(r/sigma)/r and its derivatives. The damped potential replaces the
 * bare Coulomb 1/r to eliminate the short-range divergence that causes
 * SCC convergence failures when MM atoms approach QM atoms.
 *
 * All inputs and outputs are in atomic units (Bohr, Hartree).
 */

#ifndef GRODFTB_DAMPING_H
#define GRODFTB_DAMPING_H

#ifdef __cplusplus
extern "C" {
#endif

/* -----------------------------------------------------------------------
 * Gaussian-Damped Coulomb Function
 * ThDD:T-US-043b-1.3, ThDD:T-US-043b-9.1
 *
 * f(r, sigma) = erf(r/sigma) / r
 *
 * For r/sigma < 0.01, uses Taylor expansion (Eq. T-US-043b-9.1):
 *   f = (2 / sigma / sqrt(pi)) * [1 - u^2/3 + u^4/10 - u^6/42 + u^8/216]
 * where u = r/sigma.
 *
 * Limits:
 *   r -> 0:      f -> 2 / (sigma * sqrt(pi))
 *   r >> sigma:  f -> 1/r (bare Coulomb)
 *   sigma -> 0:  f -> 1/r (for r > 0)
 *
 * @param r      Distance in Bohr (must be >= 0)
 * @param sigma  Damping width in Bohr (must be > 0)
 * @return f(r, sigma) in Bohr^-1
 * ----------------------------------------------------------------------- */
double grodftb_damped_coulomb(double r, double sigma);

/* -----------------------------------------------------------------------
 * Radial Derivative of Damped Coulomb (positive kernel)
 * ThDD:T-US-043b-2.4, ThDD:T-US-043b-9.2
 *
 * Returns g(r, sigma) = -df/dr >= 0:
 *   g = erf(r/sigma)/r^2 - 2*exp(-r^2/sigma^2) / (sigma*sqrt(pi)*r)
 *
 * For r/sigma < 0.01, uses Taylor expansion (Eq. T-US-043b-9.4):
 *   g = (4r) / (3*sigma^3*sqrt(pi)) * [1 - 3u^2/5 + 3u^4/14 - ...]
 *
 * Limits:
 *   r -> 0:      g -> 0 (vanishes linearly)
 *   r >> sigma:  g -> 1/r^2 (bare Coulomb)
 *
 * @param r      Distance in Bohr (must be >= 0)
 * @param sigma  Damping width in Bohr (must be > 0)
 * @return g(r, sigma) in Bohr^-2, always >= 0
 * ----------------------------------------------------------------------- */
double grodftb_damped_coulomb_deriv(double r, double sigma);

/* -----------------------------------------------------------------------
 * Combined Damped Coulomb Function and Gradient Kernel
 * ThDD:T-US-043b-1.11, ThDD:T-US-043b-2.4, ThDD:T-US-043b-9.1, ThDD:T-US-043b-9.2
 *
 * Computes both the damped potential f(r) and the gradient kernel K(r):
 *   f(r) = erf(r/sigma) / r
 *   K(r) = -f'(r) / r = g(r) / r
 *
 * K(r) is the kernel used in the force inner loop:
 *   F_alpha = Q * K * (R_alpha - R'_alpha)
 *
 * This combined call avoids the caller computing K = g/r separately,
 * which would be 0/0 at r=0. The Taylor expansion handles both f and K
 * at small r simultaneously.
 *
 * For r/sigma < 0.01:
 *   f = (2/sigma/sqrt(pi)) * [1 - u^2/3 + u^4/10 - u^6/42]
 *   K = (4/3/sigma^3/sqrt(pi)) * [1 - 3u^2/5 + 3u^4/14 - u^6/6]
 *
 * @param r      Distance in Bohr (must be >= 0)
 * @param sigma  Damping width in Bohr (must be > 0)
 * @param f_out  Output: damped potential value f(r) in Bohr^-1
 * @param K_out  Output: gradient kernel K(r) = g(r)/r in Bohr^-3
 * ----------------------------------------------------------------------- */
void grodftb_damped_coulomb_and_kernel(double r, double sigma,
                                        double *f_out, double *K_out);

#ifdef __cplusplus
}
#endif

#endif /* GRODFTB_DAMPING_H */
