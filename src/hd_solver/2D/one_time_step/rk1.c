#include "one_time_step.h"

void rk1(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, double dt)
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
    //double dz = fg_prev->dz;
    //double dx = fg_prev->dx;

    // Burde jeg ikke egentlig gjøre dette til slutt? Og kanskje før man starter å løse
    solve_elliptic_equation(bg, fg_prev); // Getting p1

    // Setting updated values of p1 to fg
    for (int i = 0; i < nz_full; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            fg->p1[i][j] = fg_prev->p1[i][j];
        }
    }

    extrapolate_2D_array(fg->p1, nz_full, nz_ghost, nx);

    // First finding dt
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            //dt = C / (fg_prev->vx[i][j]/dx + fg_prev->vz[i][j]/dz);
        }
    }

    // Solving diff eqs
    for (int i = nz_ghost+1; i < nz_full - nz_ghost-1; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            fg->s1[i][j] = fg_prev->s1[i][j] + dt * rhs_ds1_dt(bg, fg_prev, i, j);
            fg->vx[i][j] = fg_prev->vx[i][j] + dt * rhs_dvx_dt(bg, fg_prev, i, j);
            fg->vz[i][j] = fg_prev->vz[i][j] + dt * rhs_dvz_dt(bg, fg_prev, i, j);
        }
    }
    
    // Solving boundaries
    for (int j = 0; j < nx; j++)
    {
        // Top boundary
        fg->s1[nz_ghost][j] = 0.0;
        fg->vx[nz_ghost][j] = rhs_dvx_dt_vertical_boundary(bg, fg_prev, nz_ghost, j);
        fg->vz[nz_ghost][j] = 0.0;

        // Bottom boundary
        fg->s1[nz_full-nz_ghost-1][j] = 0.0;
        fg->vx[nz_full-nz_ghost-1][j] = rhs_dvx_dt_vertical_boundary(bg, fg_prev, nz_full-nz_ghost-1, j);
        fg->vz[nz_full-nz_ghost-1][j] = 0.0;
    }
    

    // Extrapolate
    extrapolate_2D_array(fg->s1, nz_full, nz_ghost, nx);
    extrapolate_2D_array(fg->vx, nz_full, nz_ghost, nx);
    extrapolate_2D_array(fg->vz, nz_full, nz_ghost, nx);


    // Solving algebraic equations. Eq of state should be in separate function
    double c_p;
    double total_density, total_temperature, total_pressure;
    for (int i = 0; i < nz_full; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            total_density = bg->rho0[i] + fg->rho1[i][j];
            total_temperature = bg->T0[i] + fg->T1[i][j];
            total_pressure = bg->p0[i] + fg->p1[i][j];
            c_p = total_pressure/(1-1/GAMMA) / (total_density*total_temperature);

            fg->T1[i][j] = (fg->s1[i][j]/c_p + (GAMMA-1.0)/GAMMA * fg->p1[i][j]/bg->p0[i]) * bg->T0[i];
            
            fg->rho1[i][j] = (fg->p1[i][j]/bg->p0[i] - fg->T1[i][j]/bg->T0[i]) * bg->rho0[i];
        }
    }
}