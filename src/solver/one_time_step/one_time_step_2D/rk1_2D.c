#include "one_time_step_2D.h"
#include "solver/equation_of_state/equation_of_state_2D/equation_of_state_2D.h"
#include "solver/first_law_of_thermodynamics/first_law_of_thermodynamics_2D/first_law_of_thermodynamics_2D.h"
#include "solver/boundary/boundary_2D/boundary_2D.h"
#include "array_utilities/array_memory_management/array_memory_management.h"

#include "spacial_derivatives_module/derivatives_2D/derivatives_2D.h"
#include <math.h>
#include "global_constants.h"
#include "global_parameters.h"

static inline int periodic_boundary(int i, int limit) {
    return (i + limit-1) % (limit-1);}

FLOAT_P rk1_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables2D *precalc, FLOAT_P dt_last, bool first_timestep)
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
    FLOAT_P dt = get_dt_2D(fg_prev, bg, grid_info, dt_last, first_timestep);

    FLOAT_P reduced_dt;
    // Picking smallest dt from all processes
    MPI_Allreduce(&dt, &reduced_dt, 1, MPI_FLOAT_P, MPI_MIN, MPI_COMM_WORLD);
    dt = reduced_dt;

    FLOAT_P rhs_s1, rhs_vy, rhs_vz;
    // Solving diff eqs for entire grid and boundary

    #if REMOVE_AVG_VZ_X == 1
        calculate_vz_horizontal_average_2D(fg_prev, grid_info, precalc);
    #endif
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {   
        for (int j = 0; j < ny; j++)
        {
            rhs_s1 = rhs_ds1_dt_2D(bg, fg_prev, grid_info, precalc, i, j);
            rhs_vy = rhs_dvy_dt_2D(bg, fg_prev, grid_info, precalc, i, j);
            rhs_vz = rhs_dvz_dt_2D(bg, fg_prev, grid_info, precalc, i, j);
            calculate_viscosity_thermal_diffusivity_coeffs_2D(bg, fg_prev, grid_info, precalc, i, j);
            step_entropy_velocity_2D(fg, fg_prev, grid_info, precalc, rhs_s1, rhs_vy, rhs_vz, dt, i, j);
        }
    }

    apply_vertical_boundary_damping_2D(fg, bg, grid_info, mpi_info, precalc, dt);

    update_vertical_boundary_entropy_velocity_2D(fg, grid_info, mpi_info);

    // Need these for T1 and rho1
    for (int i = 0; i < nz_full; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->p1[i][j] = fg_prev->p1[i][j];
        }
    }

    // Solving algebraic equations.
    first_law_of_thermodynamics_2D(fg, bg, grid_info);
    equation_of_state_2D(fg, bg, grid_info, mpi_info);
    
    // Solving elliptic equation
    solve_elliptic_equation_2D(bg, fg_prev, fg, grid_info, mpi_info, precalc); // Getting p1
    update_vertical_boundary_pressure_2D(fg, bg, grid_info, mpi_info);
    
    return dt;
}