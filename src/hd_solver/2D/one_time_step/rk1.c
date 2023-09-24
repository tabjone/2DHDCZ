#include "one_time_step.h"
#include <float.h>

double rk1(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, bool first_run)
{
    // Steps for each iteration:
    // Solve elliptic equation
    // Extrapolate p1, v1, s1
    // Solve diffeqs
    // Solve algebraic equations
    // Solve elliptic equation again

    int nx = grid_info->nx;
    int nz_ghost = grid_info->nz_ghost;
    int nz_full = grid_info->nz_full;
    double dz = grid_info->dz;
    double dx = grid_info->dx;

    // Finding dt
    double dt = DBL_MAX;  // Set to the maximum representable finite floating-point value
    #if UPWIND_ORDER == 1
    // First finding dt
    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            double new_dt = CFL_CUT *1.0 / (fabs(fg_prev->vx[i][j])/dx + fabs(fg_prev->vz[i][j])/dz);
            if (new_dt < dt) 
            {
                dt = new_dt;  // update dt if the new value is smaller
            }
        }
    }
    #elif UPWIND_ORDER == 2
    // This is not stable anyway
    dt = 1.0;
    #endif

    if (dt > MAX_DT)
    {
        dt = MAX_DT;
    }
    if (first_run && dt > 1.0)
    {
        dt = 1.0;
    }

    // Solving diff eqs
    for (int i = nz_ghost+1; i < nz_full - nz_ghost-1; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            fg->s1[i][j] = fg_prev->s1[i][j] + dt * rhs_ds1_dt(bg, fg_prev, grid_info, i, j);
            fg->vx[i][j] = fg_prev->vx[i][j] + dt * rhs_dvx_dt(bg, fg_prev, grid_info, i, j);
            fg->vz[i][j] = fg_prev->vz[i][j] + dt * rhs_dvz_dt(bg, fg_prev, grid_info, i, j);
        }
    }
    
    // Solving boundaries
    for (int j = 0; j < nx; j++)
    {
        // Top boundary
        fg->s1[nz_ghost][j] = 0.0;
        fg->vx[nz_ghost][j] = rhs_dvx_dt_horizontal_boundary(bg, fg_prev, grid_info, nz_ghost, j);
        fg->vz[nz_ghost][j] = 0.0;

        // Bottom boundary
        fg->s1[nz_full-nz_ghost-1][j] = 0.0;
        fg->vx[nz_full-nz_ghost-1][j] = rhs_dvx_dt_horizontal_boundary(bg, fg_prev, grid_info, nz_full-nz_ghost-1, j);
        fg->vz[nz_full-nz_ghost-1][j] = 0.0;
    }
    // Extrapolate
    extrapolate_2D_array(fg->s1, nz_full, nz_ghost, nx);
    extrapolate_2D_array(fg->vx, nz_full, nz_ghost, nx);
    extrapolate_2D_array(fg->vz, nz_full, nz_ghost, nx);

    // Burde jeg ikke egentlig gjøre dette til slutt? Og kanskje før man starter å løse
    solve_elliptic_equation(bg, fg_prev, fg, grid_info); // Getting p1
    extrapolate_2D_array(fg->p1, nz_full, nz_ghost, nx);


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

    return dt;
}