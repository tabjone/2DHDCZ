#include "one_time_step.h"

FLOAT_P rk1(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, FLOAT_P dt_last, bool first_timestep)
{
    int nx = grid_info->nx;
    int nz_ghost = grid_info->nz_ghost;
    int nz_full = grid_info->nz_full;

    FLOAT_P dt = get_dt(fg_prev, grid_info, dt_last, first_timestep);

    // Solving diff eqs
    for (int i = nz_ghost+1; i < nz_full - nz_ghost-1; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            fg->s1[i][j] = fg_prev->s1[i][j] + dt * rhs_ds1_dt(bg, fg_prev, grid_info, i, j);
            fg->vx[i][j] = fg_prev->vx[i][j] + dt * rhs_dvx_dt(bg, fg_prev, grid_info, i, j);
            fg->vz[i][j] = fg_prev->vz[i][j] + dt * rhs_dvz_dt(bg, fg_prev, grid_info, i, j);
        }
    }
    
    // Solving boundaries
    for (int j = 0; j < nx; j++)
    {
        // Bottom boundary
        fg->s1[nz_ghost][j] = 0.0;
        fg->vx[nz_ghost][j] = rhs_dvx_dt_vertical_boundary(bg, fg_prev, grid_info, nz_ghost, j);
        fg->vz[nz_ghost][j] = 0.0;

        // Top boundary
        fg->s1[nz_full-nz_ghost-1][j] = 0.0;
        fg->vx[nz_full-nz_ghost-1][j] = rhs_dvx_dt_vertical_boundary(bg, fg_prev, grid_info, nz_full-nz_ghost-1, j);
        fg->vz[nz_full-nz_ghost-1][j] = 0.0;
    }
    // Extrapolate
    extrapolate_2D_array(fg->s1, nz_full, nz_ghost, nx);
    extrapolate_2D_array(fg->vx, nz_full, nz_ghost, nx);
    extrapolate_2D_array(fg->vz, nz_full, nz_ghost, nx);

    // Solving elliptic equation
    solve_elliptic_equation(bg, fg_prev, fg, grid_info); // Getting p1
    extrapolate_2D_array(fg->p1, nz_full, nz_ghost, nx);

    // Solving algebraic equations.
    first_law_thermodynamics(fg, bg, grid_info);
    equation_of_state(fg, bg, grid_info);

    return dt;
}