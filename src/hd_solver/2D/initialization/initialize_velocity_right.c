#include "initialization.h"

void initialize_velocity_right(struct ForegroundVariables *fg, struct GridInfo *grid_info)
{
    int nz_ghost = grid_info->nz_ghost;
    int ny = grid_info->ny;
    int nz_full = grid_info->nz_full;

    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->vy[i][j] = 1000.0;
            fg->vz[i][j] = 0.0;
            fg->rho1[i][j] = 0.0;
            fg->p1[i][j] = 0.0;
            fg->T1[i][j] = 0.0;
            fg->s1[i][j] = 0.0;
        }
    }
}