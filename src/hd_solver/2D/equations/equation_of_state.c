#include "equations.h"

void equation_of_state(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo *grid_info)
{
    int nx = grid_info->nx;
    int nz_full = grid_info->nz_full;
    
    for (int i = 0; i < nz_full; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            fg->rho1[i][j] = (fg->p1[i][j]/bg->p0[i] - fg->T1[i][j]/bg->T0[i]) * bg->rho0[i];
        }
    }
}