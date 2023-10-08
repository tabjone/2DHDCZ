#include "rhs_functions.h"

#if DIMENSIONS == 2
FLOAT_P rhs_dvy_dt_2D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j)
{
    /*
    Calculates the right hands side of the momentum equation for the y-direction in 2D.

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
    FLOAT_P dz = grid_info->dz;
    
    // Creating pointers to foreground arrays
    FLOAT_P **p1 = fg->p1;
    FLOAT_P **vy = fg->vy;
    FLOAT_P **vz = fg->vz;

    // Creating pointers to background arrays
    FLOAT_P *one_over_rho0 = bg->one_over_rho0;

    // Periodic boundary conditions
    int j_minus = periodic_boundary(j-1, ny);
    int j_plus = periodic_boundary(j+1, ny);
    #if UPWIND_ORDER > 1
        int j_minus2 = periodic_boundary(j-2, ny);
        int j_plus2 = periodic_boundary(j+2, ny);
    #endif // UPWIND_ORDER

    // Calculate the derivatives
    FLOAT_P dp1_dy, dvy_dy, dvy_dz;
    #if UPWIND_ORDER == 1
        if (vy[i][j] >= 0)
        {
            dvy_dy = backward_first_derivative_first_order(vy[i][j], vy[i][j_minus], dy);
        }
        else
        {
            dvy_dy = forward_first_derivative_first_order(vy[i][j], vy[i][j_plus], dy);
        }
        if (vz[i][j] >= 0)
        {
            dvy_dz = backward_first_derivative_first_order(vy[i][j], vy[i-1][j], dz);
        }
        else
        {
            dvy_dz = forward_first_derivative_first_order(vy[i][j], vy[i+1][j], dz);
        }
    #elif UPWIND_ORDER == 2
        if (vy[i][j] >= 0.0)
        {
            dvy_dy = backward_first_derivative_second_order(vy[i][j], vy[i][j_minus], vy[i][j_minus2], dy);
        }
        else
        {
            dvy_dy = forward_first_derivative_second_order(vy[i][j], vy[i][j_plus], vy[i][j_plus2], dy);
        }
        if (vz[i][j] >= 0)
        {
            dvy_dz = backward_first_derivative_second_order(vy[i][j], vy[i-1][j], vy[i-2][j], dz);
        }
        else
        {
            dvy_dz = forward_first_derivative_second_order(vy[i][j], vy[i+1][j], vy[i+2][j], dz);
        }
    #endif // UPWIND_ORDER

    #if CENTRAL_ORDER == 2
        dp1_dy = central_first_derivative_second_order(p1[i][j_minus], p1[i][j_plus], dy);
    #endif // CENTRAL_ORDER

    #if GAS_PRESSURE_ON == 1
        rhs -= one_over_rho0[i]* dp1_dy;
    #endif // GAS_PRESSURE_ON

    // Advective term
    #if ADVECTION_ON == 1
        rhs -= vy[i][j]*dvy_dy + vz[i][j]*dvy_dz;
    #endif // ADVECTION_ON

    return rhs;
}
#endif // DIMENSIONS == 2