/*
 * SDD:specs.md:§6.2 — In-place unit conversion functions for GRO-DFTB
 *
 * All functions multiply each element of the array by the appropriate
 * conversion factor. Array length is 3*n (3 components per atom).
 */

#include "grodftb/units.h"

/* SDD:specs.md:§6.2 */
void grodftb_coords_nm_to_bohr(int n, double *coords)
{
    const int len = 3 * n;
    for (int i = 0; i < len; i++)
        coords[i] *= GRODFTB_NM_TO_BOHR;
}

/* SDD:specs.md:§6.2 */
void grodftb_coords_bohr_to_nm(int n, double *coords)
{
    const int len = 3 * n;
    for (int i = 0; i < len; i++)
        coords[i] *= GRODFTB_BOHR_TO_NM;
}

/* SDD:specs.md:§6.2 */
void grodftb_forces_au_to_gmx(int n, double *forces)
{
    const int len = 3 * n;
    for (int i = 0; i < len; i++)
        forces[i] *= GRODFTB_FORCE_AU_TO_GMX;
}

/* SDD:specs.md:§6.2 */
void grodftb_forces_gmx_to_au(int n, double *forces)
{
    const int len = 3 * n;
    for (int i = 0; i < len; i++)
        forces[i] *= GRODFTB_FORCE_GMX_TO_AU;
}
