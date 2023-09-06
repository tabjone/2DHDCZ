#include "extrapolation.h"

void extrapolate_2D(struct ForegroundVariables2D *foreground_variables)
{
    int nz_full = foreground_variables->nz_full;
    int nz_ghost = foreground_variables->nz_ghost;
    int nx = foreground_variables->nx;

    // Constant extrapolation
    #if EXTAPOLATE_GHOST_CELLS == 0
    for (int i = 0; i < nz_ghost; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            foreground_variables->p1[i][j] = foreground_variables->p1[nz_ghost][j];
            foreground_variables->rho1[i][j] = foreground_variables->rho1[nz_ghost][j];
            foreground_variables->T1[i][j] = foreground_variables->T1[nz_ghost][j];
            foreground_variables->s1[i][j] = foreground_variables->s1[nz_ghost][j];
            foreground_variables->vx[i][j] = foreground_variables->vx[nz_ghost][j];
            foreground_variables->vz[i][j] = foreground_variables->vz[nz_ghost][j];
            
            foreground_variables->p1[nz_full - nz_ghost + i][j] = foreground_variables->p1[nz_full - nz_ghost - 1][j];
            foreground_variables->rho1[nz_full - nz_ghost + i][j] = foreground_variables->rho1[nz_full - nz_ghost - 1][j];
            foreground_variables->T1[nz_full - nz_ghost + i][j] = foreground_variables->T1[nz_full - nz_ghost - 1][j];
            foreground_variables->s1[nz_full - nz_ghost + i][j] = foreground_variables->s1[nz_full - nz_ghost - 1][j];
            foreground_variables->vx[nz_full - nz_ghost + i][j] = foreground_variables->vx[nz_full - nz_ghost - 1][j];
            foreground_variables->vz[nz_full - nz_ghost + i][j] = foreground_variables->vz[nz_full - nz_ghost - 1][j];
        }
    }
    #endif
}