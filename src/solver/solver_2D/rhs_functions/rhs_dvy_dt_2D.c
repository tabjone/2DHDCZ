#include "../../spacial_derivatives_module/spacial_derivatives_module.h"
#include "global_float_precision.h"
#include "../../data_structures/background_data/background_variables_struct.h"
#include "../../data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "../../data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"

FLOAT_P rhs_dvy_dt_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct PrecalculatedVariables *precalc, int i, int j)
{
    /*
    Calculates the right hands side of the momentum equation for the y-direction in 2D.

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
        The right hand side of the momentum equation for the y-direction.
    */

    FLOAT_P rhs = 0.0; // This is the return value
    
    // Creating pointers to foreground arrays
    FLOAT_P **p1 = fg->p1;
    FLOAT_P **vy = fg->vy;
    FLOAT_P **vz = fg->vz;

    // Creating pointers to background arrays
    FLOAT_P *one_over_rho0 = precalc->one_over_rho0;

    // Calculate the derivatives
    FLOAT_P dp1_dy = central_first_derivative_y_2D(p1, i, j, grid_info->ny, precalc->one_over_2dy);
    FLOAT_P dvy_dy = upwind_first_derivative_y_2D(vy, vy, i, j, grid_info->ny, precalc->one_over_dy, precalc->one_over_2dy);
    FLOAT_P dvy_dz = upwind_first_derivative_z_2D(vy, vz, i, j, precalc->one_over_dz, precalc->one_over_2dz);
 
    #if GAS_PRESSURE_ON == 1
        rhs -= one_over_rho0[i]* dp1_dy;
    #endif // GAS_PRESSURE_ON

    // Advective term
    #if ADVECTION_ON == 1
        rhs -= vy[i][j]*dvy_dy + vz[i][j]*dvy_dz;
    #endif // ADVECTION_ON

    #if VISCOSITY_ON == 1

        FLOAT_P dd_vy_dz = central_second_derivative_z_2D(vy, i, j, precalc->one_over_dzdz);
        FLOAT_P dd_vz_dydz = central_second_derivative_yz_2D(vz, i, j, grid_info->ny, precalc->one_over_4dydz);

        rhs += precalc->VIS_COEFF_over_rho0[i]*(dd_vy_dz + dd_vz_dydz);
    #endif // VISCOSITY_ON

    return rhs;
}