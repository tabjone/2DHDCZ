#include "rhs_functions.h"

FLOAT_P rhs_dvx_dt_3D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j, int k)
{
    /*
    Calculates the right hands side of the momentum equation for the x-direction in 3D.

    Parameters
    ----------
    bg : struct
        A pointer to the BackgroundVariables struct.
    fg : struct
        A pointer to the ForegroundVariables struct.
    grid_info : struct
        A pointer to the GridInfo struct.
    i : int
        The index of the current cell in the z-direction.
    j : int
        The index of the current cell in the y-direction.
    k : int
        The index of the current cell in the x-direction.

    Returns
    -------
    rhs : FLOAT_P
        The right hand side of the momentum equation for the x-direction.
    */

    FLOAT_P rhs = 0.0; // This is the return value

    // Getting the grid info
    int nx = grid_info->nx;
    int ny = grid_info->ny;
    FLOAT_P dx = grid_info->dx;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P dz = grid_info->dz;
    
    // Creating pointers to foreground arrays
    FLOAT_P ***rho1 = fg->rho1;
    FLOAT_P ***p1 = fg->p1;
    FLOAT_P ***vx = fg->vx;
    FLOAT_P ***vy = fg->vy;
    FLOAT_P ***vz = fg->vz;

    // Creating pointers to background arrays
    FLOAT_P *one_over_rho0 = bg->one_over_rho0;

    // Periodic boundary conditions
    int j_minus = periodic_boundary(j-1, ny);
    int j_plus = periodic_boundary(j+1, ny);
    int k_minus = periodic_boundary(k-1, nx);
    int k_plus = periodic_boundary(k+1, nx);
    #if UPWIND_ORDER > 1
        int j_minus2 = periodic_boundary(j-2, ny);
        int j_plus2 = periodic_boundary(j+2, ny);
        int k_minus2 = periodic_boundary(k-2, nx);
        int k_plus2 = periodic_boundary(k+2, nx);
    #endif // UPWIND_ORDER

    // Calculate the derivatives
    FLOAT_P dp1_dx, dvx_dx, dvx_dy, dvx_dz;
    #if UPWIND_ORDER == 1
        if (vx[i][j][k] >= 0)
        {
            dvx_dx = backward_first_derivative_first_order(vx[i][j][k], vx[i][j][k_minus], dx);
        }
        else
        {
            dvx_dx = forward_first_derivative_first_order(vx[i][j][k], vx[i][j][k_plus], dx);
        }
        if (vy[i][j][k] >= 0)
        {
            dvx_dy = backward_first_derivative_first_order(vx[i][j][k], vx[i][j_minus][k], dy);
        }
        else
        {
            dvx_dy = forward_first_derivative_first_order(vx[i][j][k], vx[i][j_plus][k], dy);
        }
        if (vz[i][j][k] >= 0)
        {
            dvx_dz = backward_first_derivative_first_order(vx[i][j][k], vx[i-1][j][k], dz);
        }
        else
        {
            dvx_dz = forward_first_derivative_first_order(vx[i][j][k], vx[i+1][j][k], dz);
        }
    #elif UPWIND_ORDER == 2
        if (vx[i][j][k] >= 0.0)
        {
            dvx_dx = backward_first_derivative_second_order(vx[i][j][k], vx[i][j][k_minus], vx[i][k][k_minus2], dx);
        }
        else
        {
            dvx_dx = forward_first_derivative_second_order(vx[i][j], vx[i][j][k_plus], vx[i][j][k_plus2], dx);
        }
        if (vy[i][j][k] >= 0.0)
        {
            dvx_dy = backward_first_derivative_second_order(vx[i][j][k], vx[i][j_minus][k], vx[i][j_minus2][k], dy);
        }
        else
        {
            dvx_dy = forward_first_derivative_second_order(vx[i][j][k], vx[i][j_plus][k], vx[i][j_plus2][k], dy);
        }
        if (vz[i][j][k] >= 0)
        {
            dvx_dz = backward_first_derivative_second_order(vx[i][j][k], vx[i-1][j][k], vx[i-2][j][k], dz);
        }
        else
        {
            dvx_dz = forward_first_derivative_second_order(vx[i][j][k], vx[i+1][j][k], vx[i+2][j][k], dz);
        }
    #endif // UPWIND_ORDER

    #if CENTRAL_ORDER == 2
        dp1_dx = central_first_derivative_second_order(p1[i][j][k_minus], p1[i][j][k_plus], dx);
    #endif // CENTRAL_ORDER

    #if GAS_PRESSURE_ON == 1
        rhs -= one_over_rho0[i]* dp1_dx;
    #endif // GAS_PRESSURE_ON

    // Advective term
    #if ADVECTION_ON == 1
        rhs -= vx[i][j][k]*dvx_dx + vy[i][j][k]*dvx_dy + vz[i][j][k]*dvx_dz;
    #endif // ADVECTION_ON

    return rhs;
}