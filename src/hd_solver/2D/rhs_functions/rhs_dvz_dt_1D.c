#include "rhs_functions.h"

FLOAT_P rhs_dvz_dt_1D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i)
{
    /*
    Calculates the right hands side of the momentum equation for the z-direction in 1D.

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
    
    Returns
    -------
    rhs : FLOAT_P
        The right hand side of the momentum equation for the z-direction.
    */

    FLOAT_P rhs = 0.0; // This is the return value

    // Getting the grid info
    FLOAT_P dz = grid_info->dz;
    
    // Creating pointers to foreground arrays
    FLOAT_P *rho1 = fg->rho1;
    FLOAT_P *p1 = fg->p1;
    FLOAT_P *vz = fg->vz;

    // Creating pointers to background arrays
    FLOAT_P *one_over_rho0 = bg->one_over_rho0;
    FLOAT_P *g = bg->g;

    // Calculate the derivatives
    FLOAT_P dp1_dz, dvz_dz;
    #if UPWIND_ORDER == 1
        if (vz[i] >= 0)
        {
            dvz_dz = backward_first_derivative_first_order(vz[i], vz[i-1], dz);
        }
        else
        {
            dvz_dz = forward_first_derivative_first_order(vz[i], vz[i+1], dz);
        }
    #elif UPWIND_ORDER == 2
        if (vz[i] >= 0)
        {
            dvz_dz = backward_first_derivative_second_order(vz[i], vz[i-1], vz[i-2], dz);
        }
        else
        {
            dvz_dz = forward_first_derivative_second_order(vz[i], vz[i+1], vz[i+2], dz);
        }
    #endif // UPWIND_ORDER

    #if CENTRAL_ORDER == 2
        dp1_dz = central_first_derivative_second_order(p1[i-1], p1[i+1], dz);
    #endif // CENTRAL_ORDER

    #if GAS_PRESSURE_ON == 1
        rhs -= one_over_rho0[i] * dp1_dz;
    #endif // GAS_PRESSURE_ON

    #if GRAVITY_ON == 1
        rhs -= one_over_rho0[i]*rho1[i] * g[i];
    #endif // GRAVITY_ON

    // Advective term
    #if ADVECTION_ON == 1
        rhs -= vz[i]*dvz_dz;
    #endif // ADVECTION_ON

    return rhs;
}