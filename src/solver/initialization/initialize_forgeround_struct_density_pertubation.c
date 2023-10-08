#include "initialization.h"

void initialize_foreground_struct_density_pertubation(struct ForegroundVariables *fg, struct BackgroundVariables *bg, struct GridInfo *grid_info)
{
    int nz = grid_info->nz;
    int nz_ghost = grid_info->nz_ghost;
    int ny = grid_info->ny;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P dz = grid_info->dz;
    int i, j;


    for (i = nz_ghost; i < nz+nz_ghost; i++)
    {
        for (j = 0; j < ny; j++)
        {
            
            fg->rho1[i][j] = gaussian((i-nz_ghost)*dz, j*dy, 0.5*dz*nz, 0.25*dy*ny, 0.1*dz*nz, 0.1*dy*ny, 1.0e-5);
            fg->p1[i][j] = GAMMA*bg->p0[i]*fg->rho1[i][j]/bg->rho0[i];
            fg->T1[i][j] = bg->T0[i]*((GAMMA-1)/GAMMA * fg->p1[i][j]/bg->p0[i]);
            fg->s1[i][j] = 0.0;
            fg->vy[i][j] = 0.0;
            fg->vz[i][j] = 0.0;
        }
    }

    extrapolate_2D(fg, grid_info);
}