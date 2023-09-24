#include "initialization.h"

void initialize_velocity_right(struct ForegroundVariables2D *fg, struct GridInfo *grid_info)
{
    int nz_ghost = grid_info->nz_ghost;
    int nx = grid_info->nx;
    int nz_full = grid_info->nz_full;

    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            fg->vx[i][j] = 1000.0;
            fg->vz[i][j] = 0.0;
            fg->rho1[i][j] = 0.0;
            fg->p1[i][j] = 0.0;
            fg->T1[i][j] = 0.0;
            fg->s1[i][j] = 0.0;
        }
    }
}