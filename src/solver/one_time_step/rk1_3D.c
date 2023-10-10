#include "one_time_step.h"

#if DIMENSIONS == 3
FLOAT_P rk1_3D(struct BackgroundVariables *bg, struct ForegroundVariables *fg_prev, struct ForegroundVariables *fg, struct GridInfo *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt_last, bool first_timestep)
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

    // Solving elliptic equation
    solve_elliptic_equation(bg, fg_prev, fg, grid_info); // Getting p1

    // Extrapolating p1 to ghost cells
    extrapolate_3D_array_constant_down(fg->p1, grid_info); // Extrapolating p1 to ghost cells below
    extrapolate_3D_array_constant_up(fg->p1, grid_info); // Extrapolating p1 to ghost cells above

    // Getting grid info
    int nx = grid_info->nx;
    int ny = grid_info->ny;
    int nz_ghost = grid_info->nz_ghost;
    int nz_full = grid_info->nz_full;

    // Calculating damping factor
    FLOAT_P damping_factor[nz_full];
    calculate_damping(damping_factor, grid_info, mpi_info);

    // Calculating dt
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
            for (int k = 0; k < nx; k++)
            {
                fg->s1[i][j][k] = fg_prev->s1[i][j][k] + dt * damping_factor[i]*rhs_ds1_dt_3D(bg, fg_prev, grid_info, i, j, k);
                fg->vx[i][j][k] = fg_prev->vx[i][j][k] + dt * rhs_dvx_dt_3D(bg, fg_prev, grid_info, i, j, k);
                fg->vy[i][j][k] = fg_prev->vy[i][j][k] + dt * rhs_dvy_dt_3D(bg, fg_prev, grid_info, i, j, k);
                fg->vz[i][j][k] = fg_prev->vz[i][j][k] + dt * damping_factor[i]*rhs_dvz_dt_3D(bg, fg_prev, grid_info, i, j, k);
            }
        }
    }

    // Updating vy for the vertical boundaries
    for (int j = 0; j < ny; j++)
    {   
        for (int k = 0; k < nx; k++)
        {
            // Bottom boundary
            if (mpi_info->has_neighbor_below)
            {
                fg->vy[nz_ghost][j][k] = rhs_dvy_dt_3D(bg, fg_prev, grid_info, nz_ghost, j, k);
                fg->vx[nz_ghost][j][k] = rhs_dvx_dt_3D(bg, fg_prev, grid_info, nz_ghost, j, k);
            }
            else
            {
                fg->vy[nz_ghost][j][k] = rhs_dvy_dt_3D_vertical_boundary(bg, fg_prev, grid_info, nz_ghost, j, k);
                fg->vx[nz_ghost][j][k] = rhs_dvx_dt_3D_vertical_boundary(bg, fg_prev, grid_info, nz_ghost, j, k);
            }

            // Top boundary
            if (mpi_info->has_neighbor_above)
            {
                fg->vy[nz_full-nz_ghost-1][j][k] = rhs_dvy_dt_3D(bg, fg_prev, grid_info, nz_full-nz_ghost-1, j, k);
                fg->vx[nz_full-nz_ghost-1][j][k] = rhs_dvx_dt_3D(bg, fg_prev, grid_info, nz_full-nz_ghost-1, j, k);
            }
            else
            {
                fg->vy[nz_full-nz_ghost-1][j][k] = rhs_dvy_dt_3D_vertical_boundary(bg, fg_prev, grid_info, nz_full-nz_ghost-1, j, k);
                fg->vx[nz_full-nz_ghost-1][j][k] = rhs_dvx_dt_3D_vertical_boundary(bg, fg_prev, grid_info, nz_full-nz_ghost-1, j, k);
            }
        }
    }
    
    // Extrapolatating s1, vy and vz to ghost cells
    extrapolate_3D_array_constant_down(fg->s1, grid_info); // Extrapolating s1 to ghost cells below
    extrapolate_3D_array_constant_up(fg->s1, grid_info); // Extrapolating s1 to ghost cells above
    extrapolate_3D_array_constant_down(fg->vy, grid_info); // Extrapolating vy to ghost cells below
    extrapolate_3D_array_constant_up(fg->vy, grid_info); // Extrapolating vy to ghost cells above
    extrapolate_3D_array_constant_down(fg->vz, grid_info); // Extrapolating vz to ghost cells below
    extrapolate_3D_array_constant_up(fg->vz, grid_info); // Extrapolating vz to ghost cells above

    // Solving algebraic equations.
    first_law_thermodynamics(fg, bg, grid_info);
    equation_of_state(fg, bg, grid_info);

    return dt;
}
#endif // DIMENSIONS == 3