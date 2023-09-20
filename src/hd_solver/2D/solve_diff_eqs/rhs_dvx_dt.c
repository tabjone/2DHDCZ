#include "solve_diff_eqs.h"

double rhs_dvx_dt(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, int i, int j)
{
    int nx = fg->nx;

    double **p1 = fg->p1;
    double **vx = fg->vx;
    double **vz = fg->vz;
    double *rho0 = bg->rho0;

    double dp1_dx, dvx_dx, dvx_dz;

    double dx = fg->dx;
    double dz = fg->dz;

    // Periodic boundary conditions
    int j_minus = periodic_boundary(j-1, nx);
    int j_minus2 = periodic_boundary(j-2, nx);
    int j_plus = periodic_boundary(j+1, nx);
    int j_plus2 = periodic_boundary(j+2, nx);

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