#include "solve_diff_eqs.h"

FLOAT_P rhs_dvx_dt(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, int i, int j)
{
    int nx = grid_info->nx;

    FLOAT_P **p1 = fg->p1;
    FLOAT_P **vx = fg->vx;
    FLOAT_P **vz = fg->vz;
    FLOAT_P *rho0 = bg->rho0;

    FLOAT_P dp1_dx, dvx_dx, dvx_dz;

    FLOAT_P dx = grid_info->dx;
    FLOAT_P dz = grid_info->dz;

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
        dvx_dx = backward_first_derivative_first_order(vx[i][j], vx[i][j_minus], dx);
    }
    else
    {
        dvx_dx = forward_first_derivative_first_order(vx[i][j], vx[i][j_plus], dx);
    }
    if (vz[i][j] >= 0)
    {
        dvx_dz = backward_first_derivative_first_order(vx[i][j], vx[i-1][j], dz);
    }
    else
    {
        dvx_dz = forward_first_derivative_first_order(vx[i][j], vx[i+1][j], dz);
    }
    #elif UPWIND_ORDER == 2
    if (vx[i][j] >= 0.0)
    {
        dvx_dx = backward_first_derivative_second_order(vx[i][j], vx[i][j_minus], vx[i][j_minus2], dx);
    }
    else
    {
        dvx_dx = forward_first_derivative_second_order(vx[i][j], vx[i][j_plus], vx[i][j_plus2], dx);
    }
    if (vz[i][j] >= 0.0)
    {
        dvx_dz = backward_first_derivative_second_order(vx[i][j], vx[i-1][j], vx[i-2][j], dz);
    }
    else
    {
        dvx_dz = forward_first_derivative_second_order(vx[i][j], vx[i+1][j], vx[i+2][j], dz);
    }
    #else
    // Print error message
    printf("Error: UPWIND_ORDER must be 1 or 2\n");
    #endif

    #if CENTRAL_ORDER == 2
    dp1_dx = central_first_derivative_second_order(p1[i][j_minus], p1[i][j_plus], dx);
    #else
    // Print error message
    printf("Error: CENTRAL_ORDER must be 2\n");
    #endif

    return -1.0/rho0[i] * dp1_dx - vx[i][j]*dvx_dx - vz[i][j]*dvx_dz;
}