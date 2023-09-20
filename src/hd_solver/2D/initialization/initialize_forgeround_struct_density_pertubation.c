#include "initialization.h"

void initialize_foreground_struct_pertubation(struct ForegroundVariables2D *fg)
{
    int nz = fg->nz;
    int nz_ghost = fg->nz_ghost;
    int nx = fg->nx;
    double dx = fg->dx;
    double dz = fg->dz;
    int i, j;


    for (i = nz_ghost; i < nz+nz_ghost; i++)
    {
        for (j = 0; j < nx; j++)
        {
            fg->p1[i][j] = 0.0;
            fg->rho1[i][j] = gaussian((i-nz_ghost)*dz, j*dx, 0.5*dz*nz, 0.5*dx*nx, 0.1*dz*nz, 0.1*dx*nx, 100.0);
            fg->T1[i][j] = 0.0;
            fg->s1[i][j] = 0.0;
            fg->vx[i][j] = 0.0;
            fg->vz[i][j] = 0.0;
        }
    }

    extrapolate_2D(fg);
}