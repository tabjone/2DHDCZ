#include "global_float_precision.h"
#include "global_parameters.h"
#include "spacial_derivatives_module/derivatives_1D/derivatives_1D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_1D/foreground_variables_struct_1D.h"
#include "data_structures/grid_info/grid_info_1D/grid_info_struct_1D.h"
#include "data_structures/precalculated_data/precalculated_data_1D/precalculated_data_struct_1D.h"

FLOAT_P rhs_ds1_dt_1D(struct BackgroundVariables *bg, struct ForegroundVariables1D *fg, struct GridInfo1D *grid_info, struct PrecalculatedVariables1D *precalc, int i)
{
    /*
    Calculates the right hands side of the entropy equation in 1D.

    Parameters
    ----------
    bg : struct
        A pointer to the BackgroundVariables struct.
    fg : struct
        A pointer to the ForegroundVariables1D struct.
    grid_info : struct
        A pointer to the GridInfo1D struct.
    precalc : struct
        A pointer to the PrecalculatedVariables struct.
    i : int
        The index of the current cell in the z-direction.
    
    Returns
    -------
    rhs : FLOAT_P
        The right hand side of the entropy equation.
    */

    FLOAT_P rhs = 0.0; // This is the return value

    // Creating pointers to foreground arrays
    FLOAT_P *vz = fg->vz;
    FLOAT_P *s1 = fg->s1;

    // Creating pointers to background arrays
    FLOAT_P *grad_s0 = bg->grad_s0;

    FLOAT_P ds1_dz = upwind_first_derivative_z_1D(s1, vz, i, precalc->one_over_dz, precalc->one_over_2dz);

    #if ADVECTION_ON == 1
        rhs -= vz[i]*(ds1_dz + grad_s0[i]);
    #endif // ADVECTION_ON

    #if THERMAL_DIFFUSIVITY_ON == 1
        FLOAT_P *grad_T0 = precalc->grad_T0;
        FLOAT_P dd_s1_ddz = central_second_derivative_z_1D(s1, i, precalc->one_over_dzdz);

        rhs += precalc->THERM_COEFF_over_T0_rho0[i] * (grad_T0[i]*ds1_dz + bg->T0[i]*dd_s1_ddz);
    #endif // THERMAL_DIFFUSIVITY_ON
    
    return rhs;
}