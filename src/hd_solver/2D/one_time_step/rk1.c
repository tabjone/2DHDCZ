#include "one_time_step.h"

void rk1(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg)
{
    double dt = 0.1;

    int nx = fg->nx;
    int nz_ghost = fg->nz_ghost;
    int nz_full = fg->nz_full;
    //int nz = fg->nz;

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
    /*
    
    */
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

    solve_elliptic_equation(bg, fg); // Getting p1

    //solve_first_law_of_thermodynamics(bg, fg); // Getting T1

    //solve_equation_of_state(bg, fg); // Getting rho1
}