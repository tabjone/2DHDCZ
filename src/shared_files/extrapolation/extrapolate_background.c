#include "extrapolation.h"
#include <stdio.h>

void extrapolate_background(struct BackgroundVariables *bg)
{
    int nz_full = bg->nz_full;
    int nz_ghost = bg->nz_ghost;

    // Constant extrapolation
    #if EXTAPOLATE_GHOST_CELLS == 0
    for (int i = 0; i < nz_ghost; i++)
    {
        bg->r[i] = bg->r[bg->nz_ghost]-(bg->nz_ghost-i)*dz*R_SUN;
        bg->p0[i] = bg->p0[nz_ghost];
        bg->rho0[i] = bg->rho0[nz_ghost];
        bg->T0[i] = bg->T0[nz_ghost];
        bg->grad_s0[i] = bg->grad_s0[nz_ghost];
        bg->g[i] = bg->g[nz_ghost];

        bg->r[bg->nz+bg->nz_ghost+i] = bg->r[bg->nz_full-bg->nz_ghost-1]+(i+1)*dz*R_SUN;
        bg->p0[nz_full - nz_ghost + i] = bg->p0[nz_full - nz_ghost - 1];
        bg->rho0[nz_full - nz_ghost + i] = bg->rho0[nz_full - nz_ghost - 1];
        bg->T0[nz_full - nz_ghost + i] = bg->T0[nz_full - nz_ghost - 1];
        bg->grad_s0[nz_full - nz_ghost + i] = bg->grad_s0[nz_full - nz_ghost - 1];
        bg->g[nz_full - nz_ghost + i] = bg->g[nz_full - nz_ghost - 1];
    }
    #else
    // Print error message
    printf("Error: Extrapolation of ghost cells not implemented\n");
    #endif
}