#include "rhs_functions.h"

FLOAT_P rhs_ds1_dt_3D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j, int k)
{
    /*
    Calculates the right hands side of the entropy equation in 3D.

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
        The right hand side of the entropy equation.
    */

    FLOAT_P rhs = 0.0; // This is the return value

    // Getting the grid info
    int nx = grid_info->nx;
    int ny = grid_info->ny;
    FLOAT_P dx = grid_info->dx;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P dz = grid_info->dz;

    // Creating pointers to foreground arrays
    FLOAT_P ***vx = fg->vx;
    FLOAT_P ***vy = fg->vy;
    FLOAT_P ***vz = fg->vz;
    FLOAT_P ***s1 = fg->s1;

    // Creating pointers to background arrays
    FLOAT_P *grad_s0 = bg->grad_s0;

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
    #endif

    // Calculate the derivatives
    FLOAT_P ds1_dx, ds1_dy, ds1_dz;
    #if UPWIND_ORDER == 1
        if (vx[i][j][k] >= 0)
        {
            ds1_dx = backward_first_derivative_first_order(s1[i][j][k], s1[i][j][k_minus], dx);
        }
        else
        {
            ds1_dx = forward_first_derivative_first_order(s1[i][j][k], s1[i][j][k_plus], dx);
        }
        if (vy[i][j][k] >= 0)
        {
            ds1_dy = backward_first_derivative_first_order(s1[i][j], s1[i][j_minus][k], dy);
        }
        else
        {
            ds1_dy = forward_first_derivative_first_order(s1[i][j], s1[i][j_plus][k], dy);
        }
        if (vz[i][j][k] >= 0)
        {
            ds1_dz = backward_first_derivative_first_order(s1[i][j], s1[i-1][j][k], dz);
        }
        else
        {
            ds1_dz = forward_first_derivative_first_order(s1[i][j], s1[i+1][j][k], dz);
        }
    #elif UPWIND_ORDER == 2
        if (vx[i][j][k] >= 0)
        {
            ds1_dx = backward_first_derivative_second_order(s1[i][j][k], s1[i][j][k_minus], s1[i][j][k_minus2], dx);
        }
        else
        {
            ds1_dx = forward_first_derivative_second_order(s1[i][j][k], s1[i][j][k_plus], s1[i][j][k_plus2], dx);
        }
        if (vy[i][j][k] >= 0)
        {
            ds1_dy = backward_first_derivative_second_order(s1[i][j], s1[i][j_minus], s1[i][j_minus2][k], dy);
        }
        else
        {
            ds1_dy = forward_first_derivative_second_order(s1[i][j], s1[i][j_plus], s1[i][j_plus2][k], dy);
        }
        if (vz[i][j] >= 0)
        {
            ds1_dz = backward_first_derivative_second_order(s1[i][j], s1[i-1][j], s1[i-2][j][k], dz);
        }
        else
        {
            ds1_dz = forward_first_derivative_second_order(s1[i][j], s1[i+1][j], s1[i+2][j][k], dz);
        }
    #endif

    #if ADVECTION_ON == 1
        rhs -= vx[i][j][k]*ds1_dx + vy[i][j][k]*ds1_dy + vz[i][j][k]*ds1_dz;
    #endif

    return rhs;
}