#include "global_float_precision.h"
#include "global_parameters.h"
#include "spacial_derivatives_module/derivatives_1D/derivatives_1D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_1D/foreground_variables_struct_1D.h"
#include "data_structures/grid_info/grid_info_1D/grid_info_struct_1D.h"
#include "data_structures/precalculated_data/precalculated_data_1D/precalculated_data_struct_1D.h"

FLOAT_P rhs_dvz_dt_1D(struct BackgroundVariables *bg, struct ForegroundVariables1D *fg, struct GridInfo1D *grid_info, struct PrecalculatedVariables1D *precalc, int i)
{
    /*
    Calculates the right hands side of the momentum equation for the z-direction in 1D.

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
        The right hand side of the momentum equation for the z-direction.
    */

    FLOAT_P rhs = 0.0; // This is the return value

    // Creating pointers to background arrays
    FLOAT_P *one_over_rho0 = precalc->one_over_rho0;
    FLOAT_P *g = bg->g;
    
    // Creating pointers to foreground arrays
    FLOAT_P *rho1 = fg->rho1;
    FLOAT_P *p1 = fg->p1;
    FLOAT_P *vz = fg->vz;

    // Calculate the derivatives
    FLOAT_P dp1_dz = central_first_derivative_z_1D(p1, i, precalc->one_over_2dz);
    FLOAT_P dvz_dz = upwind_first_derivative_z_1D(vz, vz, i, precalc->one_over_dz, precalc->one_over_2dz);

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