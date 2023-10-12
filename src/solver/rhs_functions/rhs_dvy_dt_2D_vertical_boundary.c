#include "rhs_functions.h"

FLOAT_P rhs_dvy_dt_2D_vertical_boundary(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j)
{
    /*
    Calculates the right hands side of the momentum equation for the vertical boundary for the y-direction in 2D.

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

    Returns
    -------
    rhs : FLOAT_P
        The right hand side of the momentum equation for the y-direction.
    */

    FLOAT_P rhs = 0.0; // This is the return value

    // Getting the grid info
    int ny = grid_info->ny;
    FLOAT_P dy = grid_info->dy;
    
    // Creating pointers to foreground arrays
    FLOAT_P **p1 = fg->p1;
    FLOAT_P **vy = fg->vy;

    // Creating pointers to background arrays
    FLOAT_P *one_over_rho0 = bg->one_over_rho0;

    // Calculate the derivatives
    FLOAT_P dvy_dy = upwind_first_derivative_y(vy, vy, i, j, dy, ny);
    FLOAT_P dp1_dy = central_first_derivative_y(p1, i, j, dy, ny);

    #if GAS_PRESSURE_ON == 1
        rhs -= one_over_rho0[i]* dp1_dy;
    #endif // GAS_PRESSURE_ON

    // Advective term
    #if ADVECTION_ON == 1
        rhs -= vy[i][j]*dvy_dy;
    #endif // ADVECTION_ON

    return rhs;
}