#include "one_time_step_1D.h"
#include "solver/equation_of_state/equation_of_state_1D/equation_of_state_1D.h"
#include "solver/first_law_of_thermodynamics/first_law_of_thermodynamics_1D/first_law_of_thermodynamics_1D.h"
#include "solver/boundary/boundary_1D/boundary_1D.h"

FLOAT_P rk2_1D(struct BackgroundVariables *bg, struct ForegroundVariables1D *fg_prev, struct ForegroundVariables1D *fg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables1D *precalc, FLOAT_P dt_last, bool first_timestep)
{
    /*
    Calculates the foreground at the next timestep using the RK2 method.

    Parameters
    ----------
    bg : struct
        A pointer to the BackgroundVariables struct.
    fg_prev : struct
        A pointer to the ForegroundVariables1D struct at the previous timestep.
    fg : struct
        A pointer to the ForegroundVariables1D struct at the current timestep.
    grid_info : struct
        A pointer to the GridInfo1D struct.
    mpi_info : struct
        A pointer to the MpiInfo struct.
    dt_last : FLOAT_P
        The timestep used at the previous timestep.
    first_timestep : bool
        True if this is the first timestep, false otherwise.
    */

    // Getting grid info
    int nz_ghost = grid_info->nz_ghost;
    int nz_full = grid_info->nz_full;
    int nz = grid_info->nz;

    // Calculating dt
    FLOAT_P dt = get_dt_1D(fg_prev, grid_info, dt_last, first_timestep);

    
    FLOAT_P reduced_dt;
    // Picking smallest dt from all processes
    MPI_Allreduce(&dt, &reduced_dt, 1, MPI_FLOAT_P, MPI_MIN, MPI_COMM_WORLD);
    dt = reduced_dt;

    // Slopes
    FLOAT_P **k1_s1, **k2_s1;
    FLOAT_P **k1_vz, **k2_vz;

    // Allocate memory for k1, k2
    allocate_1D_array(&k1_s1, nz_full);
    allocate_1D_array(&k2_s1, nz_full);
    allocate_1D_array(&k1_vz, nz_full);
    allocate_1D_array(&k2_vz, nz_full);

    // Inside the grid and boundaries
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
            // Calculating k1 for the grid
            k1_s1[i] = rhs_ds1_dt_1D(bg, fg_prev, grid_info, precalc, i, j);
            k1_vy[i] = rhs_dvy_dt_1D(bg, fg_prev, grid_info, precalc, i, j);
            k1_vz[i] = rhs_dvz_dt_1D(bg, fg_prev, grid_info, precalc, i, j);
    }

    // Updating fg to hold mid-calculation variables
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
            fg->s1[i] = fg_prev->s1[i] + dt * k1_s1[i];
            fg->vy[i] = fg_prev->vy[i] + dt * k1_vy[i];
            fg->vz[i] = fg_prev->vz[i] + dt * k1_vz[i];
    }

    apply_vertical_boundary_damping_1D(fg, bg, grid_info, mpi_info, dt);
    update_vertical_boundary_entropy_velocity_1D(fg, grid_info, mpi_info);

    // Need these for T1 and rho1
    for (int i = 0; i < nz_full; i++)
    {
        fg->p1[i] = fg_prev->p1[i];
    }

    // Updating mid calculation variables
    first_law_thermodynamics_1D(fg, bg, grid_info);
    equation_of_state_1D(fg, bg, grid_info);

    // Calculating p1
    solve_elliptic_equation_1D(bg, fg, fg, grid_info, mpi_info, precalc); // Getting p1
    update_vertical_boundary_pressure_1D(fg, grid_info, mpi_info);

    // Calculating k2
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
            // Calculating k2 for the grid
        k2_s1[i] = rhs_ds1_dt_1D(bg, fg, grid_info, precalc, i);
        k2_vy[i] = rhs_dvy_dt_1D(bg, fg, grid_info, precalc, i);
        k2_vz[i] = rhs_dvz_dt_1D(bg, fg, grid_info, precalc, i);
    }

    // Updating variables for grid and boundaries
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        fg->s1[i] = fg_prev->s1[i] + dt/2.0 * (k1_s1[i] + k2_s1[i]);
        fg->vy[i] = fg_prev->vy[i] + dt/2.0 * (k1_vy[i] + k2_vy[i]);
        fg->vz[i] = fg_prev->vz[i] + dt/2.0 * (k1_vz[i] + k2_vz[i]);
    }

    apply_vertical_boundary_damping_1D(fg, bg, grid_info, mpi_info, dt);

    update_vertical_boundary_entropy_velocity_1D(fg, grid_info, mpi_info);

    // Solving algebraic equations
    first_law_thermodynamics_1D(fg, bg, grid_info);
    equation_of_state_1D(fg, bg, grid_info);

    // Calculating p1
    solve_elliptic_equation_1D(bg, fg, fg, grid_info, mpi_info, precalc);
    update_vertical_boundary_pressure_1D(fg, grid_info, mpi_info);

    // Deallocating memory
    deallocate_1D_array(k1_s1);
    deallocate_1D_array(k2_s1);
    deallocate_1D_array(k1_vz);
    deallocate_1D_array(k2_vz);

    return dt;
}