#include "rhs_functions.h"
#include "global_parameters.h"

FLOAT_P rhs_dvz_dt_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct PrecalculatedVariables *precalc, int i, int j)
{
    /*
    Calculates the right hands side of the momentum equation for the z-direction in 2D.

    Parameters
    ----------
    bg : struct
        A pointer to the BackgroundVariables struct.
    fg : struct
        A pointer to the ForegroundVariables2D struct.
    grid_info : struct
        A pointer to the GridInfo2D struct.
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

    // Creating pointers to background arrays
    FLOAT_P *one_over_rho0 = precalc->one_over_rho0;
    FLOAT_P *g = bg->g;
    
    // Creating pointers to foreground arrays
    FLOAT_P **rho1 = fg->rho1;
    FLOAT_P **p1 = fg->p1;
    FLOAT_P **vy = fg->vy;
    FLOAT_P **vz = fg->vz;

    // Calculate the derivatives
    FLOAT_P dp1_dz = central_first_derivative_z(p1, precalc, i, j);
    FLOAT_P dvz_dy = upwind_first_derivative_y(vz, vy, precalc, i, j);
    FLOAT_P dvz_dz = upwind_first_derivative_z(vz, vz, precalc, i, j);

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

    #if VISCOSITY_ON == 1

        FLOAT_P dd_vz_dy = central_second_derivative_y(vz, precalc, i, j);
        FLOAT_P dd_vy_dydz = central_second_derivative_yz(vy, precalc, i, j);

        rhs += precalc->VIS_COEFF_over_rho0[i]*(dd_vz_dy + dd_vy_dydz);
    #endif // VISCOSITY_ON

    return rhs;
}