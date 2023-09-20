



/*

void rk2(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, double dt)
{
    // Steps for each iteration:
    // Solve elliptic equation
    // Extrapolate p1, v1, s1
    // Solve diffeqs
    // Solve algebraic equations
    // Solve elliptic equation again


    int nz_ghost = fg->nz_ghost;
    int nz_full = fg->nz_full;
    int nx = fg->nx;

    solve_elliptic_equation(bg, fg_prev); // Getting p1

    // This is repeated in all methods, so make functions!!!

    // Setting updated values of p1 to fg
    for (int i = 0; i < nz_full; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            fg->p1[i][j] = fg_prev->p1[i][j];
        }
    }

    // Extrapolating p1 (constant extrapolation, move this to function later with possiblity for different extrapolation scheme)
    for (int k = 0; k < nz_ghost; k++)
    {
        for (int j = 0; j < nx; j++)
        {
            fg->p1[k][j] = fg->p1[nz_ghost][j];
            fg->p1[nz_full-k-1][j] = fg->p1[nz_full-nz_ghost-1][j];
        }
    }

    double **k1_s, **k2_s;
    double **k1_vx, **k2_vx;
    double **k1_vz, **k2_vz;

    // Allocate memory for k1, k2
    allocate_2D_array(&k1_s, nz_full, nx); // These only have to be nz long, but wait untill program runs to change it
    allocate_2D_array(&k2_s, nz_full, nx);
    allocate_2D_array(&k1_vx, nz_full, nx);
    allocate_2D_array(&k2_vx, nz_full, nx);
    allocate_2D_array(&k1_vz, nz_full, nx);
    allocate_2D_array(&k2_vz, nz_full, nx);

    // Just implementing RK4 for s1 for now to understand implementation

    double dt = 0.1;
    // Calculating k1
    for (int i = nz_ghost+1; i < nz_full - nz_ghost-1; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            k1_s[i][j] = rhs_ds1_dt(bg, fg_prev, i, j);
            k1_vx[i][j] = rhs_dvx_dt(bg, fg_prev, i, j);
            k1_vz[i][j] = rhs_dvz_dt(bg, fg_prev, i, j);

            // This is what we feed into the rhs for finding k2
            fg->s1[i][j] = fg_prev->s1[i][j] + dt * k1_s[i][j];
            fg->vx[i][j] = fg_prev->vx[i][j] + dt * k1_vx[i][j];
            fg->vz[i][j] = fg_prev->vz[i][j] + dt * k1_vz[i][j];
        }
    }
    // Calculating k2
    for (int i = nz_ghost+1; i < nz_full - nz_ghost-1; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            k2_s[i][j] = rhs_ds1_dt(bg, fg, i, j);
            k2_vx[i][j] = rhs_dvx_dt(bg, fg, i, j);
            k2_vz[i][j] = rhs_dvz_dt(bg, fg, i, j);
        }
    }

    // Find k1, k2 ect before for boundary before updating variables

    // Requirement: the vertical velocity at the top and bottom boundaries must be zero
    // Requirement: the horizontal velocity gradient at the top and bottom boundaries must be zero


    // Updating variables
    for (int i = nz_ghost+1; i < nz_full - nz_ghost-1; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            fg->s1[i][j] = fg_prev->s1[i][j] + dt/2.0 * (k1_s[i][j] + k2_s[i][j]);
            fg->vx[i][j] = fg_prev->vx[i][j] + dt/2.0 * (k1_vx[i][j] + k2_vx[i][j]);
            fg->vz[i][j] = fg_prev->vz[i][j] + dt/2.0 * (k1_vz[i][j] + k2_vz[i][j]);
        }
    }




    // Implement top and bottom boundary conditions
    // at [nz_ghost][j] and [nz_full-nz_ghost-1][j]

    //Update foreground variables
    //Euler for ex: fg = fg_prev + rhs_dvx_dt_ * dt;
    
}

}
*/