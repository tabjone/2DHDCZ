#include "global_float_precision.h"
#include "global_parameters.h"
#include "spacial_derivatives_module/derivatives_1D/derivatives_1D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_1D/foreground_variables_struct_1D.h"
#include "data_structures/grid_info/grid_info_1D/grid_info_struct_1D.h"
#include "data_structures/precalculated_data/precalculated_data_1D/precalculated_data_struct_1D.h"

FLOAT_P rhs_elliptic_eq_1D(struct BackgroundVariables *bg, struct ForegroundVariables1D *fg, struct GridInfo1D *grid_info, struct PrecalculatedVariables1D *precalc, int i)
{
    /*
    Calculates the right hand side of the elliptic equation in 1D.

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
        The right hand side of the elliptic equation.
    */

    FLOAT_P rhs = 0.0; // This is the return value

    // Creating pointers to foreground arrays
    FLOAT_P *rho1 = fg->rho1;
    FLOAT_P *vz = fg->vz;

    // Creating pointers to background arrays
    FLOAT_P *rho0 = bg->rho0;
    FLOAT_P *g = bg->g;
    FLOAT_P *grad_g = precalc->grad_g;
    FLOAT_P *grad_rho0 = precalc->grad_rho0;

    // Calculate the derivatives
    // First derivatives
    FLOAT_P dvz_dz = central_first_derivative_z_1D(vz, i, precalc->one_over_2dz);
    FLOAT_P drho1_dz = central_first_derivative_z_1D(rho1, i, precalc->one_over_2dz);

    // Second derivatives
    FLOAT_P dd_vz_ddz = central_second_derivative_z_1D(vz, i, precalc->one_over_dzdz);

    #if GRAVITY_ON == 1
    {
        rhs -= g[i]*drho1_dz + rho1[i][j]*grad_g[i];
    }
    #endif // GRAVITY_ON

    #if ADVECTION_ON == 1
    {
        rhs -= rho0[i] * (vz[i][j]*dd_vz_ddz + dvz_dz_sqrd)
               +grad_rho0[i] * vz[i][j]*dvz_dz;
    }
    #endif // ADVECTION_ON

    return rhs;
}