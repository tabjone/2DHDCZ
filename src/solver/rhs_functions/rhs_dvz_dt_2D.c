#include "rhs_functions.h"
#include "global_parameters.h"

FLOAT_P rhs_dvz_dt_2D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j)
{
    /*
    Calculates the right hands side of the momentum equation for the z-direction in 2D.

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
        The right hand side of the momentum equation for the z-direction.
    */

    FLOAT_P rhs = 0.0; // This is the return value

    // Getting the grid info
    int nz_full = grid_info->nz_full;
    int ny = grid_info->ny;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P dz = grid_info->dz;

    // Creating pointers to background arrays
    FLOAT_P *one_over_rho0 = bg->one_over_rho0;
    FLOAT_P *g = bg->g;
    
    // Creating pointers to foreground arrays
    FLOAT_P **rho1 = fg->rho1;
    FLOAT_P **p1 = fg->p1;
    FLOAT_P **vy = fg->vy;
    FLOAT_P **vz = fg->vz;

    // Calculate the derivatives
    FLOAT_P dp1_dz = central_first_derivative_z(p1, i, j, dz, nz_full);
    FLOAT_P dvz_dy = upwind_first_derivative_y(vz, vy, i, j, dy, ny);
    FLOAT_P dvz_dz = upwind_first_derivative_z(vz, vz, i, j, dz, nz_full);

    #if GAS_PRESSURE_ON == 1
        rhs -= one_over_rho0[i] * dp1_dz;
    #endif // GAS_PRESSURE_ON

    #if GRAVITY_ON == 1
        rhs -= one_over_rho0[i]*rho1[i][j] * g[i];
    #endif // GRAVITY_ON

    // Advective term
    #if ADVECTION_ON == 1
        rhs -= vy[i][j]*dvz_dy + vz[i][j]*dvz_dz;
    #endif // ADVECTION_ON

    return rhs;
}