#include "equations_3D.h"

void equation_of_state_3D(struct ForegroundVariables3D *fg, struct BackgroundVariables *bg, struct GridInfo3D *grid_info)
{
    // Getting grid info
    int nx = grid_info->nx;
    int ny = grid_info->ny;
    int nz_full = grid_info->nz_full;

    for (int i = 0; i < nz_full; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                fg->rho1[i][j][k] = (fg->p1[i][j][k]/bg->p0[i] - fg->T1[i][j][k]/bg->T0[i]) * bg->rho0[i];
            }
        }
    }
}