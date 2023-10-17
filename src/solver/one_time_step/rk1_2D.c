#include "one_time_step.h"

FLOAT_P rk1_2D(struct BackgroundVariables *bg, struct ForegroundVariables *fg_prev, struct ForegroundVariables *fg, struct GridInfo *grid_info, FLOAT_P dt_last, bool first_timestep)
{
    /*
    Calculates the foreground at the next timestep using the RK1 method.

    Parameters
    ----------
    bg : struct
        A pointer to the BackgroundVariables struct.
    fg_prev : struct
        A pointer to the ForegroundVariables struct at the previous timestep.
    fg : struct
        A pointer to the ForegroundVariables struct at the current timestep.
    grid_info : struct
        A pointer to the GridInfo struct.
    mpi_info : struct
        A pointer to the MpiInfo struct.
    dt_last : FLOAT_P
        The timestep used at the previous timestep.
    first_timestep : bool
        True if this is the first timestep, false otherwise.
    */

    // Getting grid info
    int ny = grid_info->ny;
    int nz_ghost = grid_info->nz_ghost;
    int nz_full = grid_info->nz_full;

    // Calculating damping factor
    FLOAT_P damping_factor[nz_full];
    calculate_damping(damping_factor, grid_info);

    // Calculating dt
    FLOAT_P dt = get_dt(fg_prev, grid_info, dt_last, first_timestep);

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
    #if VERTICAL_BOUNDARY_TYPE != 2
        for (int j = 0; j < ny; j++)
        {
            fg->vy[nz_ghost][j] = rhs_dvy_dt_2D_vertical_boundary(bg, fg_prev, grid_info, nz_ghost, j);
            fg->vy[nz_full-nz_ghost-1][j] = rhs_dvy_dt_2D_vertical_boundary(bg, fg_prev, grid_info, nz_full-nz_ghost-1, j);
        }
    #endif // VERTICAL_BOUNDARY_TYPE

    #if VERTICAL_BOUNDARY_TYPE == 2
        // Periodic boundary conditions
        periodic_boundary_2D(fg->vy, grid_info);
        periodic_boundary_2D(fg->vz, grid_info);
        periodic_boundary_2D(fg->s1, grid_info);

    #else
        // Extrapolatating s1, vy and vz to ghost cells
        extrapolate_2D_array_down(fg->s1, grid_info); // Extrapolating s1 to ghost cells below
        extrapolate_2D_array_up(fg->s1, grid_info); // Extrapolating s1 to ghost cells above
        extrapolate_2D_array_down(fg->vy, grid_info); // Extrapolating vy to ghost cells below
        extrapolate_2D_array_up(fg->vy, grid_info); // Extrapolating vy to ghost cells above
        extrapolate_2D_array_down(fg->vz, grid_info); // Extrapolating vz to ghost cells below
        extrapolate_2D_array_up(fg->vz, grid_info); // Extrapolating vz to ghost cells above
    #endif // VERTICAL_BOUNDARY_TYPE

    // Solving elliptic equation
    solve_elliptic_equation(bg, fg_prev, fg, grid_info); // Getting p1

    #if VERTICAL_BOUNDARY_TYPE == 2
        // Periodic boundary conditions
        periodic_boundary_2D(fg->p1, grid_info);
    // Extrapolating p1 to ghost cells
    #else
        extrapolate_2D_array_down(fg->p1, grid_info); // Extrapolating p1 to ghost cells below
        extrapolate_2D_array_up(fg->p1, grid_info); // Extrapolating p1 to ghost cells above
    #endif // VERTICAL_BOUNDARY_TYPE

    // Solving algebraic equations.
    first_law_thermodynamics(fg, bg, grid_info);
    equation_of_state(fg, bg, grid_info);

    return dt;
}