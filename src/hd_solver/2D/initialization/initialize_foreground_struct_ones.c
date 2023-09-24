#include "initialization.h"

void initialize_foreground_struct_ones(struct ForegroundVariables2D *fg, struct GridInfo *grid_info)
{
    int nz_ghost = grid_info->nz_ghost;
    int nx = grid_info->nx;
    int nz_full = grid_info->nz_full;

    for (int i = nz_ghost+1; i < nz_full - nz_ghost-1; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            fg->p1[i][j] = 1.0;
            fg->rho1[i][j] = 1.0;
            fg->T1[i][j] = 1.0;
            fg->s1[i][j] = 1.0;
            fg->vx[i][j] = 1.0;
            fg->vz[i][j] = 1.0;
        }
    }
    // boundary values should be zero
    for (int j = 0; j < nx; j++)
    {
        // Top boundary
        fg->p1[nz_ghost][j] = 0.0;
        fg->rho1[nz_ghost][j] = 0.0;
        fg->T1[nz_ghost][j] = 0.0;
        fg->s1[nz_ghost][j] = 0.0;
        fg->vx[nz_ghost][j] = 0.0;
        fg->vz[nz_ghost][j] = 0.0;

        // Bottom boundary
        fg->p1[nz_full-nz_ghost-1][j] = 0.0;
        fg->rho1[nz_full-nz_ghost-1][j] = 0.0;
        fg->T1[nz_full-nz_ghost-1][j] = 0.0;
        fg->s1[nz_full-nz_ghost-1][j] = 0.0;
        fg->vx[nz_full-nz_ghost-1][j] = 0.0;
        fg->vz[nz_full-nz_ghost-1][j] = 0.0;
    }

    extrapolate_2D(fg, grid_info);
}