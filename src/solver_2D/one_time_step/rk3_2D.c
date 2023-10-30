#include "one_time_step.h"

FLOAT_P rk3_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt_last, bool first_timestep)
{
    /*
    Calculates the foreground at the next timestep using the RK3 method.

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
    calculate_damping(damping_factor, bg, grid_info);

    // Calculating dt
    FLOAT_P dt = get_dt(fg_prev, grid_info, dt_last, first_timestep);

    // Using the fg struct to store mid-calculation variables. Filling these with fg_prev values.
    for (int i = 0; i < nz_full; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->T1[i][j] = fg_prev->T1[i][j];
            fg->rho1[i][j] = fg_prev->rho1[i][j];
        }
    }

    // Slopes
    FLOAT_P **k1_s1, **k2_s1, **k3_s1;
    FLOAT_P **k1_vy, **k2_vy, **k3_vy;
    FLOAT_P **k1_vz, **k2_vz, **k3_vz;

    // Allocate memory for k1, k2, k3
    allocate_2D_array(&k1_s1, nz_full, ny);
    allocate_2D_array(&k2_s1, nz_full, ny);
    allocate_2D_array(&k3_s1, nz_full, ny);
    allocate_2D_array(&k1_vy, nz_full, ny);
    allocate_2D_array(&k2_vy, nz_full, ny);
    allocate_2D_array(&k3_vy, nz_full, ny);
    allocate_2D_array(&k1_vz, nz_full, ny);
    allocate_2D_array(&k2_vz, nz_full, ny);
    allocate_2D_array(&k3_vz, nz_full, ny);

    // Inside the grid
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            // Calculating k1 for the grid
            k1_s1[i][j] = rhs_ds1_dt_2D(bg, fg_prev, grid_info, i, j);
            k1_vy[i][j] = rhs_dvy_dt_2D(bg, fg_prev, grid_info, i, j);
            k1_vz[i][j] = rhs_dvz_dt_2D(bg, fg_prev, grid_info, i, j);
        }
    }

    // Boundary
    for (int j = 0; j < ny; j++)
    {
        k1_vy[nz_ghost][j] = rhs_dvy_dt_2D_vertical_boundary(bg, fg_prev, grid_info, nz_ghost, j);
        k1_vy[nz_full-nz_ghost-1][j] = rhs_dvy_dt_2D_vertical_boundary(bg, fg_prev, grid_info, nz_full-nz_ghost-1, j);
        k1_vz[nz_ghost][j] = 0.0;
        k1_vz[nz_full-nz_ghost-1][j] = 0.0;
        k1_s1[nz_ghost][j] = 0.0;
        k1_s1[nz_full-nz_ghost-1][j] = 0.0;
    }

    // Updating fg to hold mid-calculation variables
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            // Using fg to store mid-calculation variables.
            fg->s1[i][j] = fg_prev->s1[i][j] + dt/2.0 * k1_s1[i][j];
            fg->vy[i][j] = fg_prev->vy[i][j] + dt/2.0 * k1_vy[i][j];
            fg->vz[i][j] = fg_prev->vz[i][j] + dt/2.0 * k1_vz[i][j];
        }
    }

    // Extrapolatating s1, vy and vz to ghost cells
    extrapolate_2D_array_down(fg->s1, nz_ghost, ny); // Extrapolating s1 to ghost cells below
    extrapolate_2D_array_up(fg->s1, nz_full, nz_ghost, ny); // Extrapolating s1 to ghost cells above
    extrapolate_2D_array_down(fg->vy, nz_ghost, ny); // Extrapolating vy to ghost cells below
    extrapolate_2D_array_up(fg->vy, nz_full, nz_ghost, ny); // Extrapolating vy to ghost cells above
    extrapolate_2D_array_down(fg->vz, nz_ghost, ny); // Extrapolating vz to ghost cells below
    extrapolate_2D_array_up(fg->vz, nz_full, nz_ghost, ny); // Extrapolating vz to ghost cells above

    // Calculating pressure
    solve_elliptic_equation(bg, fg, fg, grid_info, mpi_info);
    extrapolate_2D_array_constant_down(fg->p1, nz_ghost, ny); // Extrapolating p1 to ghost cells below
    extrapolate_2D_array_constant_up(fg->p1, nz_full, nz_ghost, ny); // Extrapolating p1 to ghost cells above

    // Updating mid-calculation variables
    first_law_thermodynamics(fg, bg, grid_info);
    equation_of_state(fg, bg, grid_info);

    // Calculating k2 inside the grid
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            k2_s1[i][j] = rhs_ds1_dt_2D(bg, fg, grid_info, i, j);
            k2_vy[i][j] = rhs_dvy_dt_2D(bg, fg, grid_info, i, j);
            k2_vz[i][j] = rhs_dvz_dt_2D(bg, fg, grid_info, i, j);
        }
    }

    // Boundary
    for (int j = 0; j < ny; j++)
    {
        k2_vy[nz_ghost][j] = rhs_dvy_dt_2D_vertical_boundary(bg, fg, grid_info, nz_ghost, j);
        k2_vy[nz_full-nz_ghost-1][j] = rhs_dvy_dt_2D_vertical_boundary(bg, fg, grid_info, nz_full-nz_ghost-1, j);
        k2_vz[nz_ghost][j] = 0.0;
        k2_vz[nz_full-nz_ghost-1][j] = 0.0;
        k2_s1[nz_ghost][j] = 0.0;
        k2_s1[nz_full-nz_ghost-1][j] = 0.0;
    }

    // Updating fg to hold mid-calculation variables
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->s1[i][j] = fg_prev->s1[i][j] - dt * k1_s1[i][j] + 2.0 * dt * k2_s1[i][j];
            fg->vy[i][j] = fg_prev->vy[i][j] - dt * k1_vy[i][j] + 2.0 * dt * k2_vy[i][j];
            fg->vz[i][j] = fg_prev->vz[i][j] - dt * k1_vz[i][j] + 2.0 * dt * k2_vz[i][j];
        }
    }

    // Extrapolatating s1, vy and vz to ghost cells
    extrapolate_2D_array_down(fg->s1, nz_ghost, ny); // Extrapolating s1 to ghost cells below
    extrapolate_2D_array_up(fg->s1, nz_full, nz_ghost, ny); // Extrapolating s1 to ghost cells above
    extrapolate_2D_array_down(fg->vy, nz_ghost, ny); // Extrapolating vy to ghost cells below
    extrapolate_2D_array_up(fg->vy, nz_full, nz_ghost, ny); // Extrapolating vy to ghost cells above
    extrapolate_2D_array_down(fg->vz, nz_ghost, ny); // Extrapolating vz to ghost cells below
    extrapolate_2D_array_up(fg->vz, nz_full, nz_ghost, ny); // Extrapolating vz to ghost cells above

    // Calculating pressure
    solve_elliptic_equation(bg, fg, fg, grid_info, mpi_info);
    extrapolate_2D_array_constant_down(fg->p1, nz_ghost, ny); // Extrapolating p1 to ghost cells below
    extrapolate_2D_array_constant_up(fg->p1, nz_full, nz_ghost, ny); // Extrapolating p1 to ghost cells above

    // Updating mid-calculation variables
    first_law_thermodynamics(fg, bg, grid_info);
    equation_of_state(fg, bg, grid_info);

    // Calculating k3 inside the grid
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            k3_s1[i][j] = rhs_ds1_dt_2D(bg, fg, grid_info, i, j);
            k3_vy[i][j] = rhs_dvy_dt_2D(bg, fg, grid_info, i, j);
            k3_vz[i][j] = rhs_dvz_dt_2D(bg, fg, grid_info, i, j);
        }
    }

    // Boundary
    for (int j = 0; j < ny; j++)
    {
        k3_vy[nz_ghost][j] = rhs_dvy_dt_2D_vertical_boundary(bg, fg, grid_info, nz_ghost, j);
        k3_vy[nz_full-nz_ghost-1][j] = rhs_dvy_dt_2D_vertical_boundary(bg, fg, grid_info, nz_full-nz_ghost-1, j);
        k3_vz[nz_ghost][j] = 0.0;
        k3_vz[nz_full-nz_ghost-1][j] = 0.0;
        k3_s1[nz_ghost][j] = 0.0;
        k3_s1[nz_full-nz_ghost-1][j] = 0.0;
    }

    // Finding mean of s1
    FLOAT_P s1_mean = 0.0;
    FLOAT_P tau = 60.0*60.0*4;
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            s1_mean += fg_prev->s1[i][j];
        }
    }
    s1_mean /= (nz_full - 2*nz_ghost)*ny;

    // Updating variables
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->s1[i][j] = damping_factor[i]*(fg_prev->s1[i][j] + dt/6.0 * (k1_s1[i][j] + 4.0*k2_s1[i][j] + k3_s1[i][j])) + s1_mean*(1-damping_factor[i])*(1-dt/tau);;
            fg->vy[i][j] = damping_factor[i]*(fg_prev->vy[i][j] + dt/6.0 * (k1_vy[i][j] + 4.0*k2_vy[i][j] + k3_vy[i][j]));
            fg->vz[i][j] = damping_factor[i]*(fg_prev->vz[i][j] + dt/6.0 * (k1_vz[i][j] + 4.0*k2_vz[i][j] + k3_vz[i][j]));
        }
    }

    // Extrapolating variables to ghost cells
    extrapolate_2D_array_down(fg->s1, nz_ghost, ny); // Extrapolating s1 to ghost cells below
    extrapolate_2D_array_up(fg->s1, nz_full, nz_ghost, ny); // Extrapolating s1 to ghost cells above
    extrapolate_2D_array_down(fg->vy, nz_ghost, ny); // Extrapolating vy to ghost cells below
    extrapolate_2D_array_up(fg->vy, nz_full, nz_ghost, ny); // Extrapolating vy to ghost cells above
    extrapolate_2D_array_down(fg->vz, nz_ghost, ny); // Extrapolating vz to ghost cells below
    extrapolate_2D_array_up(fg->vz, nz_full, nz_ghost, ny); // Extrapolating vz to ghost cells above

    // Solving elliptic equation
    solve_elliptic_equation(bg, fg, fg, grid_info, mpi_info); // Getting p1

    // Damping pressure
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->p1[i][j] = damping_factor[i]*fg->p1[i][j];
        }
    }

    // Extrapolating p1 to ghost cells
    extrapolate_2D_array_down(fg->p1, nz_ghost, ny); // Extrapolating p1 to ghost cells below
    extrapolate_2D_array_up(fg->p1, nz_full, nz_ghost, ny); // Extrapolating p1 to ghost cells above

    // Solving algebraic equations
    first_law_thermodynamics(fg, bg, grid_info);
    equation_of_state(fg, bg, grid_info);

    // Deallocating memory
    deallocate_2D_array(k1_s1);
    deallocate_2D_array(k2_s1);
    deallocate_2D_array(k3_s1);
    deallocate_2D_array(k1_vy);
    deallocate_2D_array(k2_vy);
    deallocate_2D_array(k3_vy);
    deallocate_2D_array(k1_vz);
    deallocate_2D_array(k2_vz);
    deallocate_2D_array(k3_vz);

    return dt;
}