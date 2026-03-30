/*
 * SDD:specs.md:S8.1 -- Quintic switching function for cutoff embedding
 * ThDD:T-US-039-2.1 -- S(u) = 1 - 10u^3 + 15u^4 - 6u^5
 * ThDD:T-US-039-3.1 -- S'(u) = -30u^2(1-u)^2
 *
 * The switching function smoothly reduces electrostatic interactions from
 * full strength to zero as the QM-MM distance crosses the switching region.
 * This eliminates force discontinuities at the cutoff boundary, ensuring
 * energy conservation in NVE simulations.
 *
 * Parameters:
 *   r_off = r_cut       (outer boundary, cutoff radius)
 *   r_on  = r_cut - w   (inner boundary, where switching begins)
 *   w     = switch_width (width of the switching region)
 *   u     = (r - r_on) / w (dimensionless coordinate, 0 <= u <= 1)
 *
 * Boundary conditions (from specs.md Section 8.1):
 *   S(0) = 1,  S(1) = 0      (value continuity)
 *   S'(0) = 0, S'(1) = 0     (first derivative continuity)
 *   S''(0) = 0, S''(1) = 0   (second derivative continuity)
 *
 * Units: All inputs/outputs are in atomic units (Bohr, Hartree).
 */

#ifndef GRODFTB_SWITCHING_H
#define GRODFTB_SWITCHING_H

#ifdef __cplusplus
extern "C" {
#endif

/* -----------------------------------------------------------------------
 * Switching Function S(u)
 * ThDD:T-US-039-2.1
 * S(u) = 1 - 10u^3 + 15u^4 - 6u^5
 *
 * Evaluated using factored Horner form for efficiency (T-US-039-N.1):
 * S(u) = 1 - u^3 * (6u^2 - 15u + 10)
 *
 * @param u Reduced coordinate: u = (r - r_on) / w, where w = r_off - r_on.
 *          Valid range: 0 <= u <= 1. Values outside this range are allowed
 *          but should be handled by the caller (S=1 for u<0, S=0 for u>1).
 * @return S(u) in the range [0, 1], with S(0)=1 and S(1)=0.
 * ----------------------------------------------------------------------- */
double grodftb_switch_func(double u);

/* -----------------------------------------------------------------------
 * First Derivative of Switching Function dS/du
 * ThDD:T-US-039-3.1
 * S'(u) = -30u^2(1-u)^2
 *
 * The factored form guarantees:
 *   - S'(0) = 0 (factor u^2)
 *   - S'(1) = 0 (factor (1-u)^2)
 *   - S'(u) <= 0 for all u in [0,1] (monotonically decreasing S)
 *
 * Maximum magnitude: |S'(0.5)| = 1.875 at the midpoint.
 *
 * @param u Reduced coordinate: u = (r - r_on) / w.
 * @return dS/du (dimensionless), always <= 0 for u in [0,1].
 * ----------------------------------------------------------------------- */
double grodftb_switch_deriv(double u);

/* -----------------------------------------------------------------------
 * Combined Switching Function and Derivative with Region Handling
 * ThDD:T-US-039-N.3
 *
 * Evaluates S(r) and dS/dr with explicit handling of three regions:
 *   - r <= r_on:           S = 1, dS/dr = 0 (full interaction)
 *   - r_on < r < r_off:    polynomial evaluation
 *   - r >= r_off:          S = 0, dS/dr = 0 (excluded)
 *
 * This function computes the derivative with respect to r (not u):
 *   dS/dr = (dS/du) * (du/dr) = S'(u) / w
 *
 * @param r     Distance between QM and MM atoms (Bohr).
 * @param r_on  Inner boundary of switching region: r_on = r_cut - w (Bohr).
 * @param w     Width of switching region: w = r_off - r_on (Bohr).
 *              Must be > 0. For w = 0 (hard cutoff), do not call this function.
 * @param S_out Output: switching function value S(r) in [0, 1].
 * @param dSdr_out Output: derivative dS/dr (Bohr^-1). Negative in switching region.
 * ----------------------------------------------------------------------- */
void grodftb_switch_func_and_deriv(double r, double r_on, double w,
                                    double *S_out, double *dSdr_out);

/* -----------------------------------------------------------------------
 * Second Derivative of Switching Function d^2S/du^2
 * ThDD:T-US-039-4.1
 * S''(u) = -60u(1-u)(1-2u)
 *
 * Properties:
 *   - S''(0) = 0 (factor u)
 *   - S''(1) = 0 (factor 1-u)
 *   - S''(0.5) = 0 (inflection point, factor 1-2u)
 *
 * @param u Reduced coordinate: u = (r - r_on) / w.
 * @return d^2S/du^2 (dimensionless).
 * ----------------------------------------------------------------------- */
double grodftb_switch_deriv2(double u);

#ifdef __cplusplus
}
#endif

#endif /* GRODFTB_SWITCHING_H */
