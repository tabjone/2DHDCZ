#include "equations.h"
#include "global_parameters.h"

void first_law_thermodynamics(struct ForegroundVariables *fg, struct BackgroundVariables *bg, struct GridInfo *grid_info)
{
    int ny = grid_info->ny;
    int nz_full = grid_info->nz_full;

    FLOAT_P r_star = K_B / (MU * M_U);
    FLOAT_P c_p;

    for (int i = 0; i < nz_full; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            c_p = r_star /(1.0-1.0/GAMMA);
            fg->T1[i][j] = bg->T0[i]*(fg->s1[i][j]/c_p + (GAMMA-1)/GAMMA * fg->p1[i][j]/bg->p0[i]);
        }
    }
}