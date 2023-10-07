#include "one_time_step.h"
#include "global_parameters.h"

FLOAT_P rk3(struct BackgroundVariables *bg, struct ForegroundVariables *fg_prev, struct ForegroundVariables *fg, struct GridInfo *grid_info, FLOAT_P dt_last, bool first_timestep)
{

    int nz_ghost = grid_info->nz_ghost;
    int nz_full = grid_info->nz_full;
    int ny = grid_info->ny;

    FLOAT_P dt = get_dt(fg_prev, grid_info, dt_last, first_timestep);
    #if MPI_ON == 1
        // Picking smallest dt from all processes
        MPI_Allreduce(&dt, &dt, 1, MPI_FLOAT_P, MPI_MIN, MPI_COMM_WORLD);
    #endif // MPI_ON

    FLOAT_P **k1_s, **k2_s, **k3_s;
    FLOAT_P **k1_vy, **k2_vy, **k3_vy;
    FLOAT_P **k1_vz, **k2_vz, **k3_vz;

    // Allocate memory for k1, k2, k3
    allocate_2D_array(&k1_s, nz_full, ny); // These only have to be nz long, but wait untill program runs to change it
    allocate_2D_array(&k2_s, nz_full, ny);
    allocate_2D_array(&k3_s, nz_full, ny);
    allocate_2D_array(&k1_vy, nz_full, ny);
    allocate_2D_array(&k2_vy, nz_full, ny);
    allocate_2D_array(&k3_vy, nz_full, ny);
    allocate_2D_array(&k1_vz, nz_full, ny);
    allocate_2D_array(&k2_vz, nz_full, ny);
    allocate_2D_array(&k3_vz, nz_full, ny);

    // Calculating k1
    for (int i = nz_ghost+1; i < nz_full - nz_ghost-1; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            k1_s[i][j] = rhs_ds1_dt_2D(bg, fg_prev, grid_info, i, j);
            k1_vy[i][j] = rhs_dvy_dt_2D(bg, fg_prev, grid_info, i, j);
            k1_vz[i][j] = rhs_dvz_dt_2D(bg, fg_prev, grid_info, i, j);

            // This is what we feed into the rhs for finding k2
            fg->s1[i][j] = fg_prev->s1[i][j] + dt/2.0 * k1_s[i][j];
            fg->vy[i][j] = fg_prev->vy[i][j] + dt/2.0 * k1_vy[i][j];
            fg->vz[i][j] = fg_prev->vz[i][j] + dt/2.0 * k1_vz[i][j];
        }
    }

    // Calculating k1 for the boundaries
    for (int j = 0; j < ny; j++)
    {
        // Bottom boundary
        k1_s[nz_ghost][j] = 0.0;
        k1_vy[nz_ghost][j] = 0.0;
        //rhs_dvy_dt_vertical_boundary(bg, fg_prev, grid_info, nz_ghost, j);
        k1_vz[nz_ghost][j] = 0.0;

        // Top boundary
        k1_s[nz_full-nz_ghost-1][j] = 0.0;
        k1_vy[nz_full-nz_ghost-1][j] = 0.0;
        //rhs_dvy_dt_vertical_boundary(bg, fg_prev, grid_info, nz_full-nz_ghost-1, j);
        k1_vz[nz_full-nz_ghost-1][j] = 0.0;

        // Bottom boundary
        fg->s1[nz_ghost][j] = fg_prev->s1[nz_ghost][j] + dt/2.0 * k1_s[nz_ghost][j];
        fg->vy[nz_ghost][j] = fg_prev->vy[nz_ghost][j] + dt/2.0 * k1_vy[nz_ghost][j];
        fg->vz[nz_ghost][j] = fg_prev->vz[nz_ghost][j] + dt/2.0 * k1_vz[nz_ghost][j];

        // Top boundary
        fg->s1[nz_full-nz_ghost-1][j] = fg_prev->s1[nz_full-nz_ghost-1][j] + dt/2.0 * k1_s[nz_full-nz_ghost-1][j];
        fg->vy[nz_full-nz_ghost-1][j] = fg_prev->vy[nz_full-nz_ghost-1][j] + dt/2.0 * k1_vy[nz_full-nz_ghost-1][j];
        fg->vz[nz_full-nz_ghost-1][j] = fg_prev->vz[nz_full-nz_ghost-1][j] + dt/2.0 * k1_vz[nz_full-nz_ghost-1][j];
    }

    extrapolate_2D_array(fg->s1, nz_full, nz_ghost, ny);
    extrapolate_2D_array(fg->vy, nz_full, nz_ghost, ny);
    extrapolate_2D_array(fg->vz, nz_full, nz_ghost, ny);

    // Using the fg struct to store mid-calculation variables. So need to fill these with the previous fg values for p1, T1 and rho1
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->p1[i][j] = fg_prev->p1[i][j];
            fg->T1[i][j] = fg_prev->T1[i][j];
            fg->rho1[i][j] = fg_prev->rho1[i][j];
        }
    }

    // Calculating k2
    for (int i = nz_ghost+1; i < nz_full - nz_ghost-1; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            k2_s[i][j] = rhs_ds1_dt_2D(bg, fg, grid_info, i, j);
            k2_vy[i][j] = rhs_dvy_dt_2D(bg, fg, grid_info, i, j);
            k2_vz[i][j] = rhs_dvz_dt_2D(bg, fg, grid_info, i, j);

            // This is what we feed into the rhs for finding k3
            fg->s1[i][j] = fg_prev->s1[i][j] - dt * k1_s[i][j] + 2.0 * dt * k2_s[i][j];
            fg->vy[i][j] = fg_prev->vy[i][j] - dt * k1_vy[i][j] + 2.0 * dt * k2_vy[i][j];
            fg->vz[i][j] = fg_prev->vz[i][j] - dt * k1_vz[i][j] + 2.0 * dt * k2_vz[i][j];
        }
    }

    // Boundaries
    for (int j = 0; j < ny; j++)
    {
        // Bottom boundary
        k2_s[nz_ghost][j] = 0.0;
        k2_vy[nz_ghost][j] = 0.0;
        //rhs_dvy_dt_vertical_boundary(bg, fg, grid_info, nz_ghost, j);
        k2_vz[nz_ghost][j] = 0.0;

        // Top boundary
        k2_s[nz_full-nz_ghost-1][j] = 0.0;
        k2_vy[nz_full-nz_ghost-1][j] = 0.0;
        //rhs_dvy_dt_vertical_boundary(bg, fg, grid_info, nz_full-nz_ghost-1, j);
        k2_vz[nz_full-nz_ghost-1][j] = 0.0;

        // Bottom boundary
        fg->s1[nz_ghost][j] = fg_prev->s1[nz_ghost][j] - dt * k1_s[nz_ghost][j] + 2.0 * dt * k2_s[nz_ghost][j];
        fg->vy[nz_ghost][j] = fg_prev->vy[nz_ghost][j] - dt * k1_vy[nz_ghost][j] + 2.0 * dt * k2_vy[nz_ghost][j];
        fg->vz[nz_ghost][j] = fg_prev->vz[nz_ghost][j] - dt * k1_vz[nz_ghost][j] + 2.0 * dt * k2_vz[nz_ghost][j];

        // Top boundary
        fg->s1[nz_full-nz_ghost-1][j] = fg_prev->s1[nz_full-nz_ghost-1][j] - dt * k1_s[nz_full-nz_ghost-1][j] + 2.0 * dt * k2_s[nz_full-nz_ghost-1][j];
        fg->vy[nz_full-nz_ghost-1][j] = fg_prev->vy[nz_full-nz_ghost-1][j] - dt * k1_vy[nz_full-nz_ghost-1][j] + 2.0 * dt * k2_vy[nz_full-nz_ghost-1][j];
        fg->vz[nz_full-nz_ghost-1][j] = fg_prev->vz[nz_full-nz_ghost-1][j] - dt * k1_vz[nz_full-nz_ghost-1][j] + 2.0 * dt * k2_vz[nz_full-nz_ghost-1][j];
    }

    extrapolate_2D_array(fg->s1, nz_full, nz_ghost, ny);
    extrapolate_2D_array(fg->vy, nz_full, nz_ghost, ny);
    extrapolate_2D_array(fg->vz, nz_full, nz_ghost, ny);

    // Calculating k3
    for (int i = nz_ghost+1; i < nz_full - nz_ghost-1; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            k3_s[i][j] = rhs_ds1_dt_2D(bg, fg, grid_info, i, j);
            k3_vy[i][j] = rhs_dvy_dt_2D(bg, fg, grid_info, i, j);
            k3_vz[i][j] = rhs_dvz_dt_2D(bg, fg, grid_info, i, j);
        }
    }

    for (int j = 0; j < ny; j++)
    {
        k3_s[nz_ghost][j] = 0.0;
        k3_vy[nz_ghost][j] = 0.0;
        //rhs_dvy_dt_vertical_boundary(bg, fg, grid_info, nz_ghost, j);
        k3_vz[nz_ghost][j] = 0.0;

        k3_s[nz_full-nz_ghost-1][j] = 0.0;
        k3_vy[nz_full-nz_ghost-1][j] = 0.0;
        //rhs_dvy_dt_vertical_boundary(bg, fg, grid_info, nz_full-nz_ghost-1, j);
        k3_vz[nz_full-nz_ghost-1][j] = 0.0;
    }

    //Updating variables
    for (int i = nz_ghost+1; i < nz_full - nz_ghost-1; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->s1[i][j] = fg_prev->s1[i][j] + dt/6.0 * (k1_s[i][j] + 4.0*k2_s[i][j] + k3_s[i][j]);
            fg->vy[i][j] = fg_prev->vy[i][j] + dt/6.0 * (k1_vy[i][j] + 4.0*k2_vy[i][j] + k3_vy[i][j]);
            fg->vz[i][j] = fg_prev->vz[i][j] + dt/6.0 * (k1_vz[i][j] + 4.0*k2_vz[i][j] + k3_vz[i][j]);
        }
    }

    // Extrapolate
    extrapolate_2D_array(fg->s1, nz_full, nz_ghost, ny);
    extrapolate_2D_array(fg->vy, nz_full, nz_ghost, ny);
    extrapolate_2D_array(fg->vz, nz_full, nz_ghost, ny);

    solve_elliptic_equation(bg, fg_prev, fg, grid_info); // Getting p1
    extrapolate_2D_array(fg->p1, nz_full, nz_ghost, ny);

    // Solving algebraic equations
    first_law_thermodynamics(fg, bg, grid_info);
    equation_of_state(fg, bg, grid_info);

    deallocate_2D_array(k1_s);
    deallocate_2D_array(k2_s);
    deallocate_2D_array(k3_s);
    deallocate_2D_array(k1_vy);
    deallocate_2D_array(k2_vy);
    deallocate_2D_array(k3_vy);
    deallocate_2D_array(k1_vz);
    deallocate_2D_array(k2_vz);
    deallocate_2D_array(k3_vz);

    return dt;
}