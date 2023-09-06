#include "functions.h"

void one_time_step(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg)
{

    //double rhs_dvx_dt_, rhs_dvz_dt_, rhs_ds1_dt_;
    //int nz = fg->nz;
    int nx = fg->nx;
    int nz_ghost = fg->nz_ghost;
    int nz_full = fg->nz_full;

    // Solve inside the grid
    for (int i = nz_ghost+1; i < nz_full - nz_ghost-1; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            //rhs_dvx_dt_ = rhs_dvx_dt(bg, fg_prev, i, j);
            //rhs_dvz_dt_ = rhs_dvz_dt(bg, fg_prev, i, j);
            //rhs_ds1_dt_ = rhs_ds1_dt(bg, fg_prev, i, j);
        }
    }

    // Implement top and bottom boundary conditions
    // at [nz_ghost][j] and [nz_full-nz_ghost-1][j]

    //Update foreground variables
    //Euler for ex: fg = fg_prev + rhs_dvx_dt_ * dt;
}