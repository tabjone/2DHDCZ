#include "solve_diff_eqs.h"

void solve_diff_eqs(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, struct ForegroundVariables2D *fg_prev, double dt)
{
    /*
    // DETTE ER JO Rk1 SCHEMET!!!!! 

    int nx = fg->nx;
    int nz_ghost = fg->nz_ghost;
    int nz_full = fg->nz_full;

    for (int i = nz_ghost+1; i < nz_full - nz_ghost-1; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            fg->s1[i][j] = fg_prev->s1[i][j] + dt * rhs_ds1_dt(bg, fg_prev, i, j);;
            fg->vx[i][j] = fg_prev->vx[i][j] + dt * rhs_dvx_dt(bg, fg_prev, i, j);;
            fg->vz[i][j] = fg_prev->vz[i][j] + dt * rhs_dvz_dt(bg, fg_prev, i, j);;
        }
    }
    

    // OR SHOULD s1 BE ZERO AT THE BOUNDARIES????

    for (int j = 0; j < nx; j++)
    {
        // Top boundary
        fg->s1[nz_ghost][j] = rhs_ds1_dt_vertical_boundary(bg, fg_prev, nz_ghost, j);
        fg->vx[nz_ghost][j] = rhs_dvx_dt_vertical_boundary(bg, fg_prev, nz_ghost, j);
        fg->vz[nz_full-nz_ghost-1][j] = 0.0;

        // Bottom boundary
        fg->s1[nz_full-nz_ghost-1][j] = rhs_ds1_dt_vertical_boundary(bg, fg_prev, nz_full-nz_ghost-1, j);
        fg->vx[nz_full-nz_ghost-1][j] = rhs_dvx_dt_vertical_boundary(bg, fg_prev, nz_full-nz_ghost-1, j);
        fg->vz[nz_full-nz_ghost-1][j] = 0.0;
    }
    */
    // MAYBE I SHOULD EXTRAPOLATE HERE
}