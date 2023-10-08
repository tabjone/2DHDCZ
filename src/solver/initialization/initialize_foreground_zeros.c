#include "initialization.h"

void initialize_foreground_struct_zeros(struct ForegroundVariables *fg, struct GridInfo *grid_info)
{
    int ny = grid_info->ny;
    int nz_full = grid_info->nz_full;

    for (int i = 0; i < nz_full; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->p1[i][j] = 0.0;
            fg->rho1[i][j] = 0.0;
            fg->T1[i][j] = 0.0;
            fg->s1[i][j] = 0.0;
            fg->vy[i][j] = 0.0;
            fg->vz[i][j] = 0.0;
        }
    }

    extrapolate_2D(fg, grid_info);
}