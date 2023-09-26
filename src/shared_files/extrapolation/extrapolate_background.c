#include "extrapolation.h"
#include <stdio.h>

void extrapolate_background(struct BackgroundVariables *bg, struct GridInfo *grid_info)
{
    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;
    FLOAT_P dz = grid_info->dz;

    // Constant extrapolation
    #if EXTAPOLATE_GHOST_CELLS == 0
    for (int i = 0; i < nz_ghost; i++)
    {
        bg->r[i] = bg->r[nz_ghost]-(nz_ghost-i)*dz;
        bg->p0[i] = bg->p0[nz_ghost];
        bg->rho0[i] = bg->rho0[nz_ghost];
        bg->T0[i] = bg->T0[nz_ghost];
        bg->grad_s0[i] = bg->grad_s0[nz_ghost];
        bg->g[i] = bg->g[nz_ghost];

        bg->r[nz_full - nz_ghost+i] = bg->r[nz_full-nz_ghost-1]+(i+1)*dz;
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