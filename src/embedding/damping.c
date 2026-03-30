/*
 * SDD:specs.md:§8.1 -- Gaussian-damped Coulomb functions for embedding
 * ThDD:T-US-043b-1.11 -- Damped embedding potential: erf(r/sigma)/r
 * ThDD:T-US-043b-2.4  -- Radial derivative of damped potential
 * ThDD:T-US-043b-9.1  -- Taylor expansion for small r/sigma
 *
 * Numerical methods (from docs/theory/US-043b/03_numerical_methods.md):
 * - Crossover at u = r/sigma < 0.01: Taylor expansion used
 * - At u = 0.01: direct evaluation retains 11+ digits
 * - Taylor with 4 terms at u < 0.01: relative error < 10^-16
 * - No overflow/underflow for physically meaningful sigma > 10^-154 Bohr
 */

#include "grodftb/damping.h"

#include <math.h>

/* ThDD:T-US-043b-9.3 -- Crossover criterion for Taylor expansion */
#define DAMPING_TAYLOR_CUTOFF 0.01

/*
 * ThDD:T-US-043b-1.3, ThDD:T-US-043b-9.1
 * Gaussian-damped Coulomb function: f(r, sigma) = erf(r/sigma) / r
 *
 * Small-r limit (Eq. T-US-043b-1.9): 2 / (sigma * sqrt(pi))
 * Taylor (Eq. T-US-043b-9.1):
 *   f = (2/sigma/sqrt(pi)) * [1 - u^2/3 + u^4/10 - u^6/42 + u^8/216]
 */
double grodftb_damped_coulomb(double r, double sigma)
{
    const double u = r / sigma;

    if (u < DAMPING_TAYLOR_CUTOFF) {
        /* ThDD:T-US-043b-9.1 -- Taylor expansion for small r/sigma
         * Avoids 0/0 indeterminate form in erf(u)/r.
         * 4 terms give relative error < 10^-16 for u < 0.01 */
        const double u2 = u * u;
        const double series = 1.0 - u2 * (1.0/3.0 - u2 * (1.0/10.0
                             - u2 * (1.0/42.0 - u2 / 216.0)));
        return (2.0 / (sigma * sqrt(M_PI))) * series;
    }

    /* ThDD:T-US-043b-1.3 -- Direct evaluation */
    return erf(u) / r;
}

/*
 * ThDD:T-US-043b-2.4, ThDD:T-US-043b-9.2
 * Radial derivative (positive kernel): g(r, sigma) = -df/dr
 *   = erf(r/sigma)/r^2 - 2*exp(-r^2/sigma^2) / (sigma*sqrt(pi)*r)
 *
 * Small-r Taylor (Eq. T-US-043b-9.4):
 *   g = (4r) / (3*sigma^3*sqrt(pi)) * [1 - 3u^2/5 + 3u^4/14 - u^6/18]
 */
double grodftb_damped_coulomb_deriv(double r, double sigma)
{
    const double u = r / sigma;

    if (u < DAMPING_TAYLOR_CUTOFF) {
        /* ThDD:T-US-043b-9.4 -- Taylor expansion
         * g vanishes linearly in r as r -> 0 */
        const double u2 = u * u;
        const double series = 1.0 - u2 * (3.0/5.0 - u2 * (3.0/14.0
                             - u2 / 18.0));
        return (4.0 * r) / (3.0 * sigma * sigma * sigma * sqrt(M_PI))
               * series;
    }

    /* ThDD:T-US-043b-2.4 -- Direct evaluation
     * g = erf(u)/r^2 - 2*exp(-u^2) / (sigma*sqrt(pi)*r) */
    const double inv_r = 1.0 / r;
    const double erf_u = erf(u);
    const double gauss = (2.0 / (sigma * sqrt(M_PI))) * exp(-u * u);
    return erf_u * inv_r * inv_r - gauss * inv_r;
}

/*
 * ThDD:T-US-043b-1.11, ThDD:T-US-043b-2.4, ThDD:T-US-043b-9.1, ThDD:T-US-043b-9.2
 * Combined: returns f(r) and K(r) = g(r)/r = -f'(r)/r
 *
 * K(r) is the gradient kernel used in force inner loops:
 *   grad_alpha = Q * K * (R_alpha - R'_alpha)
 *
 * Avoids the caller needing to compute g/r (which is 0/0 at r=0).
 * At small r, both f and K have finite Taylor limits:
 *   f -> 2/(sigma*sqrt(pi))
 *   K -> 4/(3*sigma^3*sqrt(pi))
 */
void grodftb_damped_coulomb_and_kernel(double r, double sigma,
                                        double *f_out, double *K_out)
{
    const double u = r / sigma;

    if (u < DAMPING_TAYLOR_CUTOFF) {
        /* ThDD:T-US-043b-9.1 -- Taylor for f */
        const double u2 = u * u;
        const double inv_sigma_sqrtpi = 2.0 / (sigma * sqrt(M_PI));
        const double f_series = 1.0 - u2 * (1.0/3.0 - u2 * (1.0/10.0
                               - u2 * (1.0/42.0 - u2 / 216.0)));
        *f_out = inv_sigma_sqrtpi * f_series;

        /* ThDD:T-US-043b-9.2 -- Taylor for K = g/r
         * K = (4/3) / (sigma^3 * sqrt(pi)) * [1 - 3u^2/5 + 3u^4/14 - u^6/18] */
        const double K_series = 1.0 - u2 * (3.0/5.0 - u2 * (3.0/14.0
                               - u2 / 18.0));
        *K_out = (4.0 / (3.0 * sigma * sigma * sigma * sqrt(M_PI)))
                 * K_series;
        return;
    }

    /* ThDD:T-US-043b-1.3 -- Direct f */
    const double inv_r = 1.0 / r;
    const double erf_u = erf(u);
    *f_out = erf_u * inv_r;

    /* ThDD:T-US-043b-2.4 -- Direct g, then K = g / r */
    const double gauss = (2.0 / (sigma * sqrt(M_PI))) * exp(-u * u);
    const double g = erf_u * inv_r * inv_r - gauss * inv_r;
    *K_out = g * inv_r;
}
