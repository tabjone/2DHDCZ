#include "./initialization.h"

void initialize_foreground_struct_zeros(struct ForegroundVariables2D *fg)
{
    int nz = fg->nz;
    int nz_ghost = fg->nz_ghost;
    int nx = fg->nx;

    int i, j;

    for (i = nz_ghost; i < nz+nz_ghost; i++)
    {
        for (j = 0; j < nx; j++)
        {
            fg->p1[i][j] = 0.0;
            fg->rho1[i][j] = 0.0;
            fg->T1[i][j] = 0.0;
            fg->s1[i][j] = 0.0;
            fg->vx[i][j] = 0.0;
            fg->vz[i][j] = 0.0;
        }
    }

    extrapolate_2D(fg);
}