#include "one_time_step.h"

FLOAT_P rk1_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt_last, bool first_timestep)
{
    /*
    Calculates the foreground at the next timestep using the RK1 method.

    Parameters
    ----------
    bg : struct
        A pointer to the BackgroundVariables struct.
    fg_prev : struct
        A pointer to the ForegroundVariables2D struct at the previous timestep.
    fg : struct
        A pointer to the ForegroundVariables2D struct at the current timestep.
    grid_info : struct
        A pointer to the GridInfo2D struct.
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
    #if MPI_ON == 0
        calculate_damping(damping_factor, bg, grid_info);
    #elif MPI_ON == 1
        calculate_damping_mpi(damping_factor, bg, grid_info, mpi_info);
    #endif // MPI_ON

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
            fg->s1[i][j] = (fg_prev->s1[i][j] + dt*rhs_ds1_dt_2D(bg, fg_prev, grid_info, i, j))* damping_factor[i];
            fg->vy[i][j] = fg_prev->vy[i][j] + dt*rhs_dvy_dt_2D(bg, fg_prev, grid_info, i, j);
            fg->vz[i][j] = (fg_prev->vz[i][j] + dt*rhs_dvz_dt_2D(bg, fg_prev, grid_info, i, j))* damping_factor[i];
        }
    }
    #if VERTICAL_BOUNDARY_TYPE != 2
        if (!mpi_info->has_neighbor_below)
        {
            for (int j = 0; j < ny; j++)
            {
                fg->vy[nz_ghost][j] = fg_prev->vy[nz_ghost][j] + dt*rhs_dvy_dt_2D_vertical_boundary(bg, fg_prev, grid_info, nz_ghost, j);
                fg->vz[nz_ghost][j] = 0.0;
                fg->s1[nz_ghost][j] = 0.0;
            }
            
        }
        if (!mpi_info->has_neighbor_above)
        {
            for (int j = 0; j < ny; j++)
            {
                fg->vy[nz_full-nz_ghost-1][j] = fg_prev->vy[nz_full-nz_ghost-1][j] + dt*rhs_dvy_dt_2D_vertical_boundary(bg, fg_prev, grid_info, nz_full-nz_ghost-1, j);
                fg->vz[nz_full-nz_ghost-1][j] = 0.0;
                fg->s1[nz_full-nz_ghost-1][j] = 0.0;
            }
        }
    #endif // VERTICAL_BOUNDARY_TYPE

    #if VERTICAL_BOUNDARY_TYPE == 2
        // Periodic boundary conditions
        periodic_boundary_2D(fg->vy, grid_info);
        periodic_boundary_2D(fg->vz, grid_info);
        periodic_boundary_2D(fg->s1, grid_info);

    #else
        // Extrapolatating s1, vy and vz to ghost cells
        extrapolate_2D_array_down(fg->s1, nz_ghost, ny); // Extrapolating s1 to ghost cells below
        extrapolate_2D_array_up(fg->s1, nz_full, nz_ghost, ny); // Extrapolating s1 to ghost cells above
        extrapolate_2D_array_down(fg->vy, nz_ghost, ny); // Extrapolating vy to ghost cells below
        extrapolate_2D_array_up(fg->vy, nz_full, nz_ghost, ny); // Extrapolating vy to ghost cells above
        extrapolate_2D_array_down(fg->vz, nz_ghost, ny); // Extrapolating vz to ghost cells below
        extrapolate_2D_array_up(fg->vz, nz_full, nz_ghost, ny); // Extrapolating vz to ghost cells above
    #endif // VERTICAL_BOUNDARY_TYPE

    #if MPI_ON == 1
        // Sending and receiving ghost cells
        communicate_2D_ghost_above_below(fg->s1, grid_info, mpi_info);
        communicate_2D_ghost_above_below(fg->vy, grid_info, mpi_info);
        communicate_2D_ghost_above_below(fg->vz, grid_info, mpi_info);
    #endif // MPI_ON
    
    // Solving elliptic equation
    solve_elliptic_equation(bg, fg_prev, fg, grid_info, mpi_info); // Getting p1
    
    #if VERTICAL_BOUNDARY_TYPE == 2
        // Periodic boundary conditions
        periodic_boundary_2D(fg->p1, grid_info);
    // Extrapolating p1 to ghost cells
    #else
        extrapolate_2D_array_down(fg->p1, nz_ghost, ny); // Extrapolating p1 to ghost cells below
        extrapolate_2D_array_up(fg->p1, nz_full, nz_ghost, ny); // Extrapolating p1 to ghost cells above
    #endif // VERTICAL_BOUNDARY_TYPE

    #if MPI_ON == 1
        // Sending and receiving ghost cells
        communicate_2D_ghost_above_below(fg->p1, grid_info, mpi_info);
    #endif // MPI_ON
    
    // Solving algebraic equations.
    first_law_thermodynamics(fg, bg, grid_info);
    equation_of_state(fg, bg, grid_info);

    return dt;
}