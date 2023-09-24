#include "initialization.h"

void initialize_foreground_struct_density_pertubation(struct ForegroundVariables2D *fg, struct GridInfo *grid_info)
{
    int nz = grid_info->nz;
    int nz_ghost = grid_info->nz_ghost;
    int nx = grid_info->nx;
    double dx = grid_info->dx;
    double dz = grid_info->dz;
    int i, j;


    for (i = nz_ghost; i < nz+nz_ghost; i++)
    {
        for (j = 0; j < nx; j++)
        {
            fg->p1[i][j] = gaussian((i-nz_ghost)*dz, j*dx, 0.5*dz*nz, 0.25*dx*nx, 0.1*dz*nz, 0.1*dx*nx, 1.0e-1);
            fg->rho1[i][j] = 0.0;
            fg->T1[i][j] = 0.0;
            fg->s1[i][j] = 0.0;
            fg->vx[i][j] = 0.0;
            fg->vz[i][j] = 0.0;
        }
    }

    extrapolate_2D(fg, grid_info);
}