#include "rhs_functions.h"

#if DIMENSIONS == 1
FLOAT_P rhs_elliptic_eq_1D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i)
{
    /*
    Calculates the right hand side of the elliptic equation in 1D.

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
        The right hand side of the elliptic equation.
    */

    FLOAT_P rhs = 0.0; // This is the return value

    // Getting the grid info
    int ny = grid_info->ny;
    FLOAT_P dz = grid_info->dz;

    // Creating pointers to foreground arrays
    FLOAT_P *rho1 = fg->rho1;
    FLOAT_P *vz = fg->vz;

    // Creating pointers to background arrays
    FLOAT_P *rho0 = bg->rho0;
    FLOAT_P *grad_g = bg->grad_g;
    FLOAT_P *g = bg->g;

    // Calculate the derivatives
    FLOAT_P dd_vz_ddz, drho1_dz;

    #if CENTRAL_ORDER == 2
        dd_vz_ddz = central_second_derivative_second_order(vz[i], vz[i-1], vz[i+1], dz);
        drho1_dz = central_first_derivative_second_order(rho1[i-1], rho1[i+1], dz);
    #endif // CENTRAL_ORDER

    #if GRAVITY_ON == 1
    {
        rhs -= g[i]*drho1_dz + rho1[i]*grad_g[i];
    }
    #endif // GRAVITY_ON

    #if ADVECTION_ON == 1
    {
        rhs -= rho0[i]*vz[i]*dd_vz_ddz;
    }
    #endif // ADVECTION_ON

    return rhs;
}
#endif // DIMENSIONS