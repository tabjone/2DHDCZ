#include "solve_diff_eqs.h"

FLOAT_P rhs_dvz_dt(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, int i, int j)
{
    int nx = grid_info->nx;
    
    FLOAT_P **rho1 = fg->rho1;
    FLOAT_P **p1 = fg->p1;
    FLOAT_P **vx = fg->vx;
    FLOAT_P **vz = fg->vz;
    FLOAT_P *rho0 = bg->rho0;
    FLOAT_P *g = bg->g;

    FLOAT_P dx = grid_info->dx;
    FLOAT_P dz = grid_info->dz;

    FLOAT_P dp1_dz, dvz_dx, dvz_dz;

    // Periodic boundary conditions
    int j_minus = periodic_boundary(j-1, nx);
    int j_plus = periodic_boundary(j+1, nx);
    #if UPWIND_ORDER > 1
    int j_minus2 = periodic_boundary(j-2, nx);
    int j_plus2 = periodic_boundary(j+2, nx);
    #endif

    #if UPWIND_ORDER == 1
    if (vx[i][j] >= 0)
    {
        dvz_dx = backward_first_derivative_first_order(vz[i][j], vz[i][j_minus], dx);
    }
    else
    {
        dvz_dx = forward_first_derivative_first_order(vz[i][j], vz[i][j_plus], dx);
    }
    if (vz[i][j] >= 0)
    {
        dvz_dz = backward_first_derivative_first_order(vz[i][j], vz[i-1][j], dz);
    }
    else
    {
        dvz_dz = forward_first_derivative_first_order(vz[i][j], vz[i+1][j], dz);
    }
    #elif UPWIND_ORDER == 2
    if (vx[i][j] >= 0.0)
    {
        dvz_dx = backward_first_derivative_second_order(vz[i][j], vz[i][j_minus], vz[i][j_minus2], dx);
    }
    else
    {
        dvz_dx = forward_first_derivative_second_order(vz[i][j], vz[i][j_plus], vz[i][j_plus2], dx);
    }

    if (vz[i][j] >= 0)
    {
        dvz_dz = backward_first_derivative_second_order(vz[i][j], vz[i-1][j], vz[i-2][j], dz);
    }
    else
    {
        dvz_dz = forward_first_derivative_second_order(vz[i][j], vz[i+1][j], vz[i+2][j], dz);
    }
    #else
    // Print error message
    printf("Error: UPWIND_ORDER must be 1 or 2\n");

    #endif

    #if CENTRAL_ORDER == 2
    dp1_dz = central_first_derivative_second_order(p1[i-1][j], p1[i+1][j], dz);
    #else
    // Print error message
    printf("Error: CENTRAL_ORDER must be 2\n");
    #endif

    return - 1.0/rho0[i] * dp1_dz - vx[i][j]*dvz_dx - vz[i][j]*dvz_dz + rho1[i][j]/rho0[i] * g[i];
}