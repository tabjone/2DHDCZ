#include "rhs_functions.h"

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
    int nz_full = grid_info->nz_full;
    int ny = grid_info->ny;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P dz = grid_info->dz;
    
    // Creating pointers to foreground arrays
    FLOAT_P **p1 = fg->p1;
    FLOAT_P **vy = fg->vy;
    FLOAT_P **vz = fg->vz;

    // Creating pointers to background arrays
    FLOAT_P *one_over_rho0 = bg->one_over_rho0;

    // Calculate the derivatives
    FLOAT_P dp1_dy = central_first_derivative_y(p1, i, j, dy, ny);
    FLOAT_P dvy_dy = upwind_first_derivative_y(vy, vy, i, j, dy, ny);
    FLOAT_P dvy_dz = upwind_first_derivative_z(vy, vz, i, j, dz, nz_full);
 
    #if GAS_PRESSURE_ON == 1
        rhs -= one_over_rho0[i]* dp1_dy;
    #endif // GAS_PRESSURE_ON

    // Advective term
    #if ADVECTION_ON == 1
        rhs -= vy[i][j]*dvy_dy + vz[i][j]*dvy_dz;
    #endif // ADVECTION_ON

    #if VISCOSITY_ON == 1
        // Periodic boundary conditions
        int j_minus = periodic_boundary(j-1, ny);
        int j_plus = periodic_boundary(j+1, ny);

        FLOAT_P dd_vy_dz = central_second_derivative_z(vy, i, j, dz, nz_full);
        FLOAT_P dd_vz_dydz = (vz[i+1][j_plus] - vz[i+1][j_minus] - vz[i-1][j_plus] + vz[i-1][j_minus])/(4.0*dy*dz);

        rhs += VISCOSITY_COEFF*one_over_rho0[i]*(dd_vy_dz + dd_vz_dydz);
    #endif // VISCOSITY_ON

    return rhs;
}