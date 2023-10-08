#include "one_time_step.h"

FLOAT_P rk1(struct BackgroundVariables *bg, struct ForegroundVariables *fg_prev, struct ForegroundVariables *fg, struct GridInfo *grid_info, FLOAT_P dt_last, bool first_timestep)
{
    int ny = grid_info->ny;
    int nz_ghost = grid_info->nz_ghost;
    int nz_full = grid_info->nz_full;
    //int nz = grid_info->nz;

    FLOAT_P damping_factor[nz_full];
    FLOAT_P damping_coeffs[5] = {0.0, 0.5, 0.75, 0.875, 0.9375};
    for (int i = nz_ghost+5; i < nz_full-nz_ghost-5+1; i++)
    {
        damping_factor[i] = 1.0;
    }
    for (int i = nz_ghost; i < nz_ghost+5; i++)
    {
        damping_factor[i] = damping_coeffs[i-nz_ghost];
        damping_factor[nz_full-i-1] = damping_coeffs[i-nz_ghost];
    }

    FLOAT_P dt = get_dt(fg_prev, grid_info, dt_last, first_timestep);
    #if MPI_ON == 1
        // Picking smallest dt from all processes
        MPI_Allreduce(&dt, &dt, 1, MPI_FLOAT_P, MPI_MIN, MPI_COMM_WORLD);
    #endif // MPI_ON

    // Solving diff eqs for entire grid and boundary
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->s1[i][j] = fg_prev->s1[i][j] + dt * damping_factor[i]*rhs_ds1_dt_2D(bg, fg_prev, grid_info, i, j);
            fg->vy[i][j] = fg_prev->vy[i][j] + dt * rhs_dvy_dt_2D(bg, fg_prev, grid_info, i, j);
            fg->vz[i][j] = fg_prev->vz[i][j] + dt * damping_factor[i]*rhs_dvz_dt_2D(bg, fg_prev, grid_info, i, j);
        }
    }
    // Updating vy for the vertical boundaries
    for (int j = 0; j < ny; j++)
    {   
        // THIS NEEDS TO BE UPDATED. JUST EXTRAPOLATE VX SO THAT ITS ANTISYMMETRICAL
        // Bottom boundary

        fg->vy[nz_ghost][j] = 0.0;
        //rhs_dvy_dt_vertical_boundary(bg, fg_prev, grid_info, nz_ghost, j);

        // Top boundary
        fg->vy[nz_full-nz_ghost-1][j] = 0.0;
        //rhs_dvy_dt_vertical_boundary(bg, fg_prev, grid_info, nz_full-nz_ghost-1, j);
    }
    // Extrapolate
    extrapolate_2D_array(fg->s1, nz_full, nz_ghost, ny);
    extrapolate_2D_array(fg->vy, nz_full, nz_ghost, ny);
    extrapolate_2D_array(fg->vz, nz_full, nz_ghost, ny);

    // Solving elliptic equation
    solve_elliptic_equation(bg, fg_prev, fg, grid_info); // Getting p1
    extrapolate_2D_array(fg->p1, nz_full, nz_ghost, ny);

    // Solving algebraic equations.
    first_law_thermodynamics(fg, bg, grid_info);
    equation_of_state(fg, bg, grid_info);

    return dt;
}