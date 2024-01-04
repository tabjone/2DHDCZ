#include "one_time_step_3D.h"
#include "solver/equation_of_state/equation_of_state_3D/equation_of_state_3D.h"
#include "solver/first_law_of_thermodynamics/first_law_of_thermodynamics_3D/first_law_of_thermodynamics_3D.h"
#include "solver/boundary/boundary_3D/boundary_3D.h"

FLOAT_P rk1_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg_prev, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables3D *precalc, FLOAT_P dt_last, bool first_timestep)
{
    /*
    Calculates the foreground at the next timestep using the RK1 method.

    Parameters
    ----------
    bg : struct
        A pointer to the BackgroundVariables struct.
    fg_prev : struct
        A pointer to the ForegroundVariables3D struct at the previous timestep.
    fg : struct
        A pointer to the ForegroundVariables3D struct at the current timestep.
    grid_info : struct
        A pointer to the GridInfo3D struct.
    mpi_info : struct
        A pointer to the MpiInfo struct.
    dt_last : FLOAT_P
        The timestep used at the previous timestep.
    first_timestep : bool
        True if this is the first timestep, false otherwise.
    */

    // Getting grid info
    int nx = grid_info->nx;
    int ny = grid_info->ny;
    int nz_ghost = grid_info->nz_ghost;
    int nz_full = grid_info->nz_full;

    // Calculating dt
    FLOAT_P dt = get_dt_3D(fg_prev, grid_info, dt_last, first_timestep);

    FLOAT_P reduced_dt;
    // Picking smallest dt from all processes
    MPI_Allreduce(&dt, &reduced_dt, 1, MPI_FLOAT_P, MPI_MIN, MPI_COMM_WORLD);
    dt = reduced_dt;
    

    // Solving diff eqs for entire grid and boundary
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                fg->s1[i][j][k] = fg_prev->s1[i][j][k] + dt*rhs_ds1_dt_3D(bg, fg_prev, grid_info, precalc, i, j, k);
                fg->vx[i][j][k] = fg_prev->vx[i][j][k] + dt*rhs_dvx_dt_3D(bg, fg_prev, grid_info, precalc, i, j, k);
                fg->vy[i][j][k] = fg_prev->vy[i][j][k] + dt*rhs_dvy_dt_3D(bg, fg_prev, grid_info, precalc, i, j, k);
                fg->vz[i][j][k] = fg_prev->vz[i][j][k] + dt*rhs_dvz_dt_3D(bg, fg_prev, grid_info, precalc, i, j, k);
            }
        }
    }

    apply_vertical_boundary_damping_3D(fg, bg, grid_info, mpi_info, dt);

    update_vertical_boundary_entropy_velocity_3D(fg, grid_info, mpi_info);

    // Need these for T1 and rho1
    for (int i = 0; i < nz_full; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                fg->p1[i][j][k] = fg_prev->p1[i][j][k];
            }
        }
    }

    // Solving algebraic equations.
    first_law_of_thermodynamics_3D(fg, bg, grid_info);
    equation_of_state_3D(fg, bg, grid_info);
    
    // Solving elliptic equation
    solve_elliptic_equation_3D(bg, fg_prev, fg, grid_info, mpi_info, precalc); // Getting p1
    update_vertical_boundary_pressure_3D(fg, grid_info, mpi_info);
    
    return dt;
}