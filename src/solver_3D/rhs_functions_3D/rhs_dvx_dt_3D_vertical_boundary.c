#include "rhs_functions_3D.h"

FLOAT_P rhs_dvx_dt_3D_vertical_boundary(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, int i, int j, int k)
{
    /*
    Calculates the right hands side of the momentum equation for the vertical boundary for the x-direction in 3D.

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
    
    // Creating pointers to foreground arrays
    FLOAT_P ***p1 = fg->p1;
    FLOAT_P ***vx = fg->vx;
    FLOAT_P ***vy = fg->vy;

    // Creating pointers to background arrays
    FLOAT_P *one_over_rho0 = bg->one_over_rho0;

    // Calculate the derivatives
    FLOAT_P dp1_dx = central_first_derivative_x_3D(p1, i, j, k, dx, nx);

    FLOAT_P dvx_dx = upwind_first_derivative_x_3D(vx, vx, i, j, k, dx, nx);
    FLOAT_P dvx_dy = upwind_first_derivative_y_3D(vx, vy, i, j, k, dy, ny);

    #if GAS_PRESSURE_ON == 1
        rhs -= one_over_rho0[i]* dp1_dx;
    #endif // GAS_PRESSURE_ON

    // Advective term
    #if ADVECTION_ON == 1
        rhs -= vx[i][j][k]*dvx_dx + vy[i][j][k]*dvx_dy;
    #endif // ADVECTION_ON

    return rhs;
}