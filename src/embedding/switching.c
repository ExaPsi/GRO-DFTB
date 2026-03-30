/*
 * SDD:specs.md:S8.1 -- Quintic switching function implementation
 * ThDD:T-US-039-2.1 -- S(u) = 1 - 10u^3 + 15u^4 - 6u^5
 * ThDD:T-US-039-3.1 -- S'(u) = -30u^2(1-u)^2
 * ThDD:T-US-039-4.1 -- S''(u) = -60u(1-u)(1-2u)
 *
 * Implementation notes:
 * - Factored Horner form for minimal operation count (T-US-039-N.1)
 * - Factored derivative form ensures exact zeros at boundaries (T-US-039-N.4)
 * - Three-region handling avoids unnecessary polynomial evaluation (T-US-039-N.3)
 *
 * Numerical stability (T-US-039-N.2):
 * - Near u=0: S(u) = 1 - 10u^3 + O(u^4), no cancellation
 * - Near u=1: Cancellation is exact by construction (algebraic identity)
 * - At u=0.5: Terms are comparable magnitude, exact result S(0.5) = 0.5
 *
 * Performance (T-US-039-N.11):
 * - grodftb_switch_func: 5 multiplications, 3 additions
 * - grodftb_switch_deriv: 4 multiplications, 1 addition
 * - Combined function avoids recomputing shared subexpressions
 */

#include "grodftb/switching.h"

/*
 * ThDD:T-US-039-N.1
 * Quintic switching function using factored Horner form.
 *
 * S(u) = 1 - 10u^3 + 15u^4 - 6u^5
 *      = 1 - u^3 * (10 - 15u + 6u^2)
 *      = 1 - u^3 * (6u^2 - 15u + 10)
 *
 * Operation count: 5 mult + 3 add (vs 6 mult + 4 add for direct evaluation)
 */
double grodftb_switch_func(double u)
{
    /*
     * ThDD:T-US-039-2.1
     * Evaluate S(u) = 1 - u^3 * (6u^2 - 15u + 10)
     */
    const double u2 = u * u;
    const double u3 = u2 * u;
    return 1.0 - u3 * (6.0 * u2 - 15.0 * u + 10.0);
}

/*
 * ThDD:T-US-039-N.4
 * First derivative using factored form.
 *
 * S'(u) = d/du [1 - 10u^3 + 15u^4 - 6u^5]
 *       = -30u^2 + 60u^3 - 30u^4
 *       = -30u^2 (1 - 2u + u^2)
 *       = -30u^2 (1-u)^2
 *
 * The factored form:
 * 1. Guarantees S'(0) = 0 (factor u^2)
 * 2. Guarantees S'(1) = 0 (factor (1-u)^2)
 * 3. Is always non-positive for u in [0,1], ensuring monotonic switching
 *
 * Operation count: 4 mult + 1 sub
 */
double grodftb_switch_deriv(double u)
{
    /*
     * ThDD:T-US-039-3.1
     * Evaluate S'(u) = -30 * u^2 * (1-u)^2
     */
    const double u2 = u * u;
    const double one_minus_u = 1.0 - u;
    const double one_minus_u_sq = one_minus_u * one_minus_u;
    return -30.0 * u2 * one_minus_u_sq;
}

/*
 * ThDD:T-US-039-4.1
 * Second derivative using factored form.
 *
 * S''(u) = d/du [-30u^2 + 60u^3 - 30u^4]
 *        = -60u + 180u^2 - 120u^3
 *        = -60u (1 - 3u + 2u^2)
 *        = -60u (1-u)(1-2u)
 *
 * Properties:
 * - S''(0) = 0 (factor u)
 * - S''(1) = 0 (factor 1-u)
 * - S''(0.5) = 0 (inflection point, factor 1-2u)
 * - S''(u) > 0 for u in (0, 0.5): concave up
 * - S''(u) < 0 for u in (0.5, 1): concave down
 *
 * Operation count: 4 mult + 2 sub
 */
double grodftb_switch_deriv2(double u)
{
    /*
     * ThDD:T-US-039-4.1
     * Evaluate S''(u) = -60 * u * (1-u) * (1-2u)
     */
    const double one_minus_u = 1.0 - u;
    const double one_minus_2u = 1.0 - 2.0 * u;
    return -60.0 * u * one_minus_u * one_minus_2u;
}

/*
 * ThDD:T-US-039-N.3
 * Combined switching function and derivative with explicit region handling.
 *
 * Three regions:
 * 1. r <= r_on: Full interaction, S = 1, dS/dr = 0
 * 2. r_on < r < r_off: Switching region, polynomial evaluation
 * 3. r >= r_off: Excluded, S = 0, dS/dr = 0
 *
 * The derivative with respect to r (not u) is:
 *   dS/dr = (dS/du) * (du/dr) = S'(u) / w
 * where w = r_off - r_on is the switching width.
 *
 * Branch prediction (T-US-039-N.3):
 * - Most MM atoms are either fully inside or outside the switching region
 * - Early exit for these cases avoids unnecessary polynomial evaluation
 */
void grodftb_switch_func_and_deriv(double r, double r_on, double w,
                                    double *S_out, double *dSdr_out)
{
    /*
     * ThDD:T-US-039-N.3
     * Region 1: r <= r_on (full interaction)
     * S = 1, dS/dr = 0
     */
    if (r <= r_on)
    {
        *S_out = 1.0;
        *dSdr_out = 0.0;
        return;
    }

    /*
     * ThDD:T-US-039-N.3
     * Region 3: r >= r_off = r_on + w (excluded)
     * S = 0, dS/dr = 0
     */
    const double r_off = r_on + w;
    if (r >= r_off)
    {
        *S_out = 0.0;
        *dSdr_out = 0.0;
        return;
    }

    /*
     * ThDD:T-US-039-N.3
     * Region 2: r_on < r < r_off (switching region)
     * Compute u, S(u), S'(u), and dS/dr = S'(u) / w
     */

    /* ThDD:T-US-039-2.1: u = (r - r_on) / w */
    const double u = (r - r_on) / w;
    const double u2 = u * u;
    const double u3 = u2 * u;
    const double one_minus_u = 1.0 - u;

    /*
     * ThDD:T-US-039-N.1
     * S(u) = 1 - u^3 * (6u^2 - 15u + 10)
     */
    *S_out = 1.0 - u3 * (6.0 * u2 - 15.0 * u + 10.0);

    /*
     * ThDD:T-US-039-N.4
     * S'(u) = -30 * u^2 * (1-u)^2
     * dS/dr = S'(u) / w
     */
    const double one_minus_u_sq = one_minus_u * one_minus_u;
    const double dSdu = -30.0 * u2 * one_minus_u_sq;
    *dSdr_out = dSdu / w;
}
