#include "one_time_step.h"

#if DIMENSIONS == 1
FLOAT_P rk3_1D(struct BackgroundVariables *bg, struct ForegroundVariables *fg_prev, struct ForegroundVariables *fg, struct GridInfo *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt_last, bool first_timestep)
{
    /*
    Calculates the foreground at the next timestep using the RK3 method.

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
    extrapolate_1D_array_down(fg->p1, grid_info); // Extrapolating p1 to ghost cells below
    extrapolate_1D_array_up(fg->p1, grid_info); // Extrapolating p1 to ghost cells above

    #if MPI_ON == 1
        // Communicate ghost cells
    #endif // MPI_ON

    // Getting grid info
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

    // Using the fg struct to store mid-calculation variables. Filling these with fg_prev values.
    for (int i = 0; i < nz_full; i++)
    {
        fg->T1[i] = fg_prev->T1[i];
        fg->rho1[i] = fg_prev->rho1[i];
    }

    // Slopes
    FLOAT_P *k1_s1, *k2_s1, *k3_s1;
    FLOAT_P *k1_vz, *k2_vz, *k3_vz;

    // Allocating memory for slopes
    allocate_1D_array(&k1_s1, nz_full);
    allocate_1D_array(&k2_s1, nz_full);
    allocate_1D_array(&k3_s1, nz_full);
    allocate_1D_array(&k1_vz, nz_full);
    allocate_1D_array(&k2_vz, nz_full);
    allocate_1D_array(&k3_vz, nz_full);

    // Calculating k1
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        k1_s1[i] = dt * damping_factor[i]*rhs_ds1_dt_1D(bg, fg_prev, grid_info, i);
        k1_vz[i] = dt * damping_factor[i]*rhs_dvz_dt_1D(bg, fg_prev, grid_info, i);
    }

    // Updating fg to store mid-calculation variables
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        fg->s1[i] = fg_prev->s1[i] + dt/2.0 * k1_s1[i];
        fg->vz[i] = fg_prev->vz[i] + dt/2.0 * k1_vz[i];
    }

    // Extrapolaing s1 and vz to ghost cells
    extrapolate_1D_array_down(fg->s1, grid_info); // Extrapolating s1 to ghost cells below
    extrapolate_1D_array_up(fg->s1, grid_info); // Extrapolating s1 to ghost cells above
    extrapolate_1D_array_down(fg->vz, grid_info); // Extrapolating vz to ghost cells below
    extrapolate_1D_array_up(fg->vz, grid_info); // Extrapolating vz to ghost cells above

    #if MPI_ON == 1
        // Communicate ghost cells
    #endif // MPI_ON

    // Calculating k2
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        k2_s1[i] = damping_factor[i]*rhs_ds1_dt_1D(bg, fg, grid_info, i);
        k2_vz[i] = damping_factor[i]*rhs_dvz_dt_1D(bg, fg, grid_info, i);
    }

    // Updating fg to store mid-calculation variables
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        fg->s1[i] = fg_prev->s1[i] - dt*k1_s1[i] + 2.0*dt*k2_s1[i];
        fg->vz[i] = fg_prev->vz[i] - dt*k1_vz[i] + 2.0*dt*k2_vz[i];
    }

    // Extrapolaing s1 and vz to ghost cells
    extrapolate_1D_array_down(fg->s1, grid_info); // Extrapolating s1 to ghost cells below
    extrapolate_1D_array_up(fg->s1, grid_info); // Extrapolating s1 to ghost cells above
    extrapolate_1D_array_down(fg->vz, grid_info); // Extrapolating vz to ghost cells below
    extrapolate_1D_array_up(fg->vz, grid_info); // Extrapolating vz to ghost cells above

    #if MPI_ON == 1
        // Communicate ghost cells
    #endif // MPI_ON

    // Calculating k3
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        k3_s1[i] = damping_factor[i]*rhs_ds1_dt_1D(bg, fg, grid_info, i);
        k3_vz[i] = damping_factor[i]*rhs_dvz_dt_1D(bg, fg, grid_info, i);
    }

    // Updating variables
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        fg->s1[i] = fg_prev->s1[i] + dt/6.0 * (k1_s1[i] + 4.0*k2_s1[i] + k3_s1[i]);
        fg->vz[i] = fg_prev->vz[i] + dt/6.0 * (k1_vz[i] + 4.0*k2_vz[i] + k3_vz[i]);
    }

    // Extrapolatating s1 and vz to ghost cells
    extrapolate_1D_array_down(fg->s1, grid_info); // Extrapolating s1 to ghost cells below
    extrapolate_1D_array_up(fg->s1, grid_info); // Extrapolating s1 to ghost cells above
    extrapolate_1D_array_down(fg->vz, grid_info); // Extrapolating vz to ghost cells below
    extrapolate_1D_array_up(fg->vz, grid_info); // Extrapolating vz to ghost cells above

    #if MPI_ON == 1
        // Communicate ghost cells
    #endif // MPI_ON

    // Solving algebraic equations.
    first_law_thermodynamics(fg, bg, grid_info);
    equation_of_state(fg, bg, grid_info);

    return dt;
}
#endif // DIMENSIONS == 1