#include "../../spacial_derivatives_module/spacial_derivatives_module.h"
#include "global_float_precision.h"
#include "../../data_structures/background_data/background_variables_struct.h"
#include "../../data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "../../data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"

FLOAT_P rhs_ds1_dt_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct PrecalculatedVariables *precalc, int i, int j)
{
    /*
    Calculates the right hands side of the entropy equation in 2D.

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
        The right hand side of the entropy equation.
    */

    FLOAT_P rhs = 0.0; // This is the return value

    // Creating pointers to foreground arrays
    FLOAT_P **vy = fg->vy;
    FLOAT_P **vz = fg->vz;
    FLOAT_P **s1 = fg->s1;

    // Creating pointers to background arrays
    FLOAT_P *grad_s0 = bg->grad_s0;

    FLOAT_P ds1_dy = upwind_first_derivative_y_2D(s1, vy, i, j, grid_info->ny, precalc->one_over_dy, precalc->one_over_2dy);
    FLOAT_P ds1_dz = upwind_first_derivative_z_2D(s1, vz, i, j, precalc->one_over_dz, precalc->one_over_2dz);

    #if ADVECTION_ON == 1
        rhs -= vy[i][j]*ds1_dy + vz[i][j]*ds1_dz + vz[i][j]*grad_s0[i];
    #endif // ADVECTION_ON

    #if VISCOSITY_ON == 1
        FLOAT_P dvy_dz = central_first_derivative_z_2D(vy, i, j, precalc->one_over_2dz);
        FLOAT_P dvz_dy = central_first_derivative_y_2D(vz, i, j, grid_info->ny, precalc->one_over_2dy);
        
        rhs += precalc->VIS_COEFF_over_T0_rho0[i] * (dvy_dz + dvz_dy)*(dvy_dz + dvz_dy);
    #endif // VISCOSITY_ON

    #if THERMAL_DIFFUSIVITY_ON == 1
        FLOAT_P *grad_T0 = precalc->grad_T0;
        
        FLOAT_P dd_s1_ddy = central_second_derivative_y_2D(s1, i, j, grid_info->ny, precalc->one_over_dydy);
        FLOAT_P dd_s1_ddz = central_second_derivative_z_2D(s1, i, j, precalc->one_over_dzdz);

        rhs += precalc->THERM_COEFF_over_T0_rho0[i] * (grad_T0[i]*(ds1_dy+ds1_dz) + bg->T0[i]*(dd_s1_ddy + dd_s1_ddz));
    #endif // THERMAL_DIFFUSIVITY_ON
    
    return rhs;
}