#include "one_time_step.h"

FLOAT_P rk1_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables *precalc, FLOAT_P dt_last, bool first_timestep)
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

    // Calculating dt
    FLOAT_P dt = get_dt(fg_prev, grid_info, dt_last, first_timestep);

    #if MPI_ON == 1
        FLOAT_P reduced_dt;
        // Picking smallest dt from all processes
        MPI_Allreduce(&dt, &reduced_dt, 1, MPI_FLOAT_P, MPI_MIN, MPI_COMM_WORLD);
        dt = reduced_dt;
    #endif // MPI_ON

    // Solving diff eqs for entire grid and boundary
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->s1[i][j] = fg_prev->s1[i][j] + dt*rhs_ds1_dt_2D(bg, fg_prev, grid_info, precalc, i, j);
            fg->vy[i][j] = fg_prev->vy[i][j] + dt*rhs_dvy_dt_2D(bg, fg_prev, grid_info, precalc, i, j);
            fg->vz[i][j] = fg_prev->vz[i][j] + dt*rhs_dvz_dt_2D(bg, fg_prev, grid_info, precalc, i, j);
        }
    }

    apply_vertical_boundary_damping(fg, bg, grid_info, mpi_info, dt);

    update_vertical_boundary_ghostcells_2D(fg->vy, grid_info, mpi_info);
    update_vertical_boundary_ghostcells_2D(fg->vz, grid_info, mpi_info);
    update_vertical_boundary_ghostcells_2D(fg->s1, grid_info, mpi_info);
    
    // Need these for T1 and rho1
    for (int i = 0; i < nz_full; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->p1[i][j] = fg_prev->p1[i][j];
        }
    }

    // Solving algebraic equations.
    first_law_thermodynamics(fg, bg, grid_info);
    equation_of_state(fg, bg, grid_info);
    
    // Solving elliptic equation
    solve_elliptic_equation(bg, fg_prev, fg, grid_info, mpi_info, precalc); // Getting p1
    update_vertical_boundary_ghostcells_2D(fg->p1, grid_info, mpi_info);
    
    return dt;
}