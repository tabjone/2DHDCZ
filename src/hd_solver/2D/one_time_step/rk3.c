#include "one_time_step.h"

FLOAT_P rk3(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, FLOAT_P dt_last, bool first_timestep)
{

    int nz_ghost = grid_info->nz_ghost;
    int nz_full = grid_info->nz_full;
    int nx = grid_info->nx;

    FLOAT_P dt = get_dt(fg_prev, grid_info, dt_last, first_timestep);

    FLOAT_P **k1_s, **k2_s, **k3_s;
    FLOAT_P **k1_vx, **k2_vx, **k3_vx;
    FLOAT_P **k1_vz, **k2_vz, **k3_vz;

    // Allocate memory for k1, k2, k3
    allocate_2D_array(&k1_s, nz_full, nx); // These only have to be nz long, but wait untill program runs to change it
    allocate_2D_array(&k2_s, nz_full, nx);
    allocate_2D_array(&k3_s, nz_full, nx);
    allocate_2D_array(&k1_vx, nz_full, nx);
    allocate_2D_array(&k2_vx, nz_full, nx);
    allocate_2D_array(&k3_vx, nz_full, nx);
    allocate_2D_array(&k1_vz, nz_full, nx);
    allocate_2D_array(&k2_vz, nz_full, nx);
    allocate_2D_array(&k3_vz, nz_full, nx);

    // Calculating k1
    for (int i = nz_ghost+1; i < nz_full - nz_ghost-1; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            k1_s[i][j] = rhs_ds1_dt(bg, fg_prev, grid_info, i, j);
            k1_vx[i][j] = rhs_dvx_dt(bg, fg_prev, grid_info, i, j);
            k1_vz[i][j] = rhs_dvz_dt(bg, fg_prev, grid_info, i, j);

            // This is what we feed into the rhs for finding k2
            fg->s1[i][j] = fg_prev->s1[i][j] + dt/2.0 * k1_s[i][j];
            fg->vx[i][j] = fg_prev->vx[i][j] + dt/2.0 * k1_vx[i][j];
            fg->vz[i][j] = fg_prev->vz[i][j] + dt/2.0 * k1_vz[i][j];
        }
    }

    // Calculating k1 for the boundaries
    for (int j = 0; j < nx; j++)
    {
        // Bottom boundary
        k1_s[nz_ghost][j] = 0.0;
        k1_vx[nz_ghost][j] = rhs_dvx_dt_vertical_boundary(bg, fg_prev, grid_info, nz_ghost, j);
        k1_vz[nz_ghost][j] = 0.0;

        // Top boundary
        k1_s[nz_full-nz_ghost-1][j] = 0.0;
        k1_vx[nz_full-nz_ghost-1][j] = rhs_dvx_dt_vertical_boundary(bg, fg_prev, grid_info, nz_full-nz_ghost-1, j);
        k1_vz[nz_full-nz_ghost-1][j] = 0.0;

        // Bottom boundary
        fg->s1[nz_ghost][j] = fg_prev->s1[nz_ghost][j] + dt/2.0 * k1_s[nz_ghost][j];
        fg->vx[nz_ghost][j] = fg_prev->vx[nz_ghost][j] + dt/2.0 * k1_vx[nz_ghost][j];
        fg->vz[nz_ghost][j] = fg_prev->vz[nz_ghost][j] + dt/2.0 * k1_vz[nz_ghost][j];

        // Top boundary
        fg->s1[nz_full-nz_ghost-1][j] = fg_prev->s1[nz_full-nz_ghost-1][j] + dt/2.0 * k1_s[nz_full-nz_ghost-1][j];
        fg->vx[nz_full-nz_ghost-1][j] = fg_prev->vx[nz_full-nz_ghost-1][j] + dt/2.0 * k1_vx[nz_full-nz_ghost-1][j];
        fg->vz[nz_full-nz_ghost-1][j] = fg_prev->vz[nz_full-nz_ghost-1][j] + dt/2.0 * k1_vz[nz_full-nz_ghost-1][j];
    }

    extrapolate_2D_array(fg->s1, nz_full, nz_ghost, nx);
    extrapolate_2D_array(fg->vx, nz_full, nz_ghost, nx);
    extrapolate_2D_array(fg->vz, nz_full, nz_ghost, nx);

    // Using the fg struct to store mid-calculation variables. So need to fill these with the previous fg values for p1, T1 and rho1
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            fg->p1[i][j] = fg_prev->p1[i][j];
            fg->T1[i][j] = fg_prev->T1[i][j];
            fg->rho1[i][j] = fg_prev->rho1[i][j];
        }
    }

    // Calculating k2
    for (int i = nz_ghost+1; i < nz_full - nz_ghost-1; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            k2_s[i][j] = rhs_ds1_dt(bg, fg, grid_info, i, j);
            k2_vx[i][j] = rhs_dvx_dt(bg, fg, grid_info, i, j);
            k2_vz[i][j] = rhs_dvz_dt(bg, fg, grid_info, i, j);

            // This is what we feed into the rhs for finding k3
            fg->s1[i][j] = fg_prev->s1[i][j] - dt * k1_s[i][j] + 2.0 * dt * k2_s[i][j];
            fg->vx[i][j] = fg_prev->vx[i][j] - dt * k1_vx[i][j] + 2.0 * dt * k2_vx[i][j];
            fg->vz[i][j] = fg_prev->vz[i][j] - dt * k1_vz[i][j] + 2.0 * dt * k2_vz[i][j];
        }
    }

    // Boundaries
    for (int j = 0; j < nx; j++)
    {
        // Bottom boundary
        k2_s[nz_ghost][j] = 0.0;
        k2_vx[nz_ghost][j] = rhs_dvx_dt_vertical_boundary(bg, fg, grid_info, nz_ghost, j);
        k2_vz[nz_ghost][j] = 0.0;

        // Top boundary
        k2_s[nz_full-nz_ghost-1][j] = 0.0;
        k2_vx[nz_full-nz_ghost-1][j] = rhs_dvx_dt_vertical_boundary(bg, fg, grid_info, nz_full-nz_ghost-1, j);
        k2_vz[nz_full-nz_ghost-1][j] = 0.0;

        // Bottom boundary
        fg->s1[nz_ghost][j] = fg_prev->s1[nz_ghost][j] - dt * k1_s[nz_ghost][j] + 2.0 * dt * k2_s[nz_ghost][j];
        fg->vx[nz_ghost][j] = fg_prev->vx[nz_ghost][j] - dt * k1_vx[nz_ghost][j] + 2.0 * dt * k2_vx[nz_ghost][j];
        fg->vz[nz_ghost][j] = fg_prev->vz[nz_ghost][j] - dt * k1_vz[nz_ghost][j] + 2.0 * dt * k2_vz[nz_ghost][j];

        // Top boundary
        fg->s1[nz_full-nz_ghost-1][j] = fg_prev->s1[nz_full-nz_ghost-1][j] - dt * k1_s[nz_full-nz_ghost-1][j] + 2.0 * dt * k2_s[nz_full-nz_ghost-1][j];
        fg->vx[nz_full-nz_ghost-1][j] = fg_prev->vx[nz_full-nz_ghost-1][j] - dt * k1_vx[nz_full-nz_ghost-1][j] + 2.0 * dt * k2_vx[nz_full-nz_ghost-1][j];
        fg->vz[nz_full-nz_ghost-1][j] = fg_prev->vz[nz_full-nz_ghost-1][j] - dt * k1_vz[nz_full-nz_ghost-1][j] + 2.0 * dt * k2_vz[nz_full-nz_ghost-1][j];
    }

    extrapolate_2D_array(fg->s1, nz_full, nz_ghost, nx);
    extrapolate_2D_array(fg->vx, nz_full, nz_ghost, nx);
    extrapolate_2D_array(fg->vz, nz_full, nz_ghost, nx);

    // Calculating k3
    for (int i = nz_ghost+1; i < nz_full - nz_ghost-1; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            k3_s[i][j] = rhs_ds1_dt(bg, fg, grid_info, i, j);
            k3_vx[i][j] = rhs_dvx_dt(bg, fg, grid_info, i, j);
            k3_vz[i][j] = rhs_dvz_dt(bg, fg, grid_info, i, j);
        }
    }

    for (int j = 0; j < nx; j++)
    {
        k3_s[nz_ghost][j] = 0.0;
        k3_vx[nz_ghost][j] = rhs_dvx_dt_vertical_boundary(bg, fg, grid_info, nz_ghost, j);
        k3_vz[nz_ghost][j] = 0.0;

        k3_s[nz_full-nz_ghost-1][j] = 0.0;
        k3_vx[nz_full-nz_ghost-1][j] = rhs_dvx_dt_vertical_boundary(bg, fg, grid_info, nz_full-nz_ghost-1, j);
        k3_vz[nz_full-nz_ghost-1][j] = 0.0;
    }

    //Updating variables
    for (int i = nz_ghost+1; i < nz_full - nz_ghost-1; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            fg->s1[i][j] = fg_prev->s1[i][j] + dt/6.0 * (k1_s[i][j] + 4.0*k2_s[i][j] + k3_s[i][j]);
            fg->vx[i][j] = fg_prev->vx[i][j] + dt/6.0 * (k1_vx[i][j] + 4.0*k2_vx[i][j] + k3_vx[i][j]);
            fg->vz[i][j] = fg_prev->vz[i][j] + dt/6.0 * (k1_vz[i][j] + 4.0*k2_vz[i][j] + k3_vz[i][j]);
        }
    }

    // Extrapolate
    extrapolate_2D_array(fg->s1, nz_full, nz_ghost, nx);
    extrapolate_2D_array(fg->vx, nz_full, nz_ghost, nx);
    extrapolate_2D_array(fg->vz, nz_full, nz_ghost, nx);

    solve_elliptic_equation(bg, fg_prev, fg, grid_info); // Getting p1
    extrapolate_2D_array(fg->p1, nz_full, nz_ghost, nx);

    // Solving algebraic equations
    first_law_thermodynamics(fg, bg, grid_info);
    equation_of_state(fg, bg, grid_info);

    deallocate_2D_array(k1_s);
    deallocate_2D_array(k2_s);
    deallocate_2D_array(k3_s);
    deallocate_2D_array(k1_vx);
    deallocate_2D_array(k2_vx);
    deallocate_2D_array(k3_vx);
    deallocate_2D_array(k1_vz);
    deallocate_2D_array(k2_vz);
    deallocate_2D_array(k3_vz);

    return dt;
}