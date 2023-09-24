#include "extrapolation.h"

void extrapolate_2D(struct ForegroundVariables2D *fg, struct GridInfo *grid_info)
{
    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;
    int nx = grid_info->nx;

    // Constant extrapolation
    #if EXTAPOLATE_GHOST_CELLS == 0
    for (int i = 0; i < nz_ghost; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            fg->p1[i][j] = fg->p1[nz_ghost][j];
            fg->rho1[i][j] = fg->rho1[nz_ghost][j];
            fg->T1[i][j] = fg->T1[nz_ghost][j];
            fg->s1[i][j] = fg->s1[nz_ghost][j];
            fg->vx[i][j] = fg->vx[nz_ghost][j];
            fg->vz[i][j] = fg->vz[nz_ghost][j];
            
            fg->p1[nz_full - nz_ghost + i][j] = fg->p1[nz_full - nz_ghost - 1][j];
            fg->rho1[nz_full - nz_ghost + i][j] = fg->rho1[nz_full - nz_ghost - 1][j];
            fg->T1[nz_full - nz_ghost + i][j] = fg->T1[nz_full - nz_ghost - 1][j];
            fg->s1[nz_full - nz_ghost + i][j] = fg->s1[nz_full - nz_ghost - 1][j];
            fg->vx[nz_full - nz_ghost + i][j] = fg->vx[nz_full - nz_ghost - 1][j];
            fg->vz[nz_full - nz_ghost + i][j] = fg->vz[nz_full - nz_ghost - 1][j];
        }
    }
    #endif
}