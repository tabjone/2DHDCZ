#include "global_float_precision.h"
#include "global_parameters.h"
#include "spacial_derivatives_module/derivatives_2D/derivatives_2D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "data_structures/precalculated_data/precalculated_data_2D/precalculated_data_struct_2D.h"

FLOAT_P rhs_elliptic_eq_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct PrecalculatedVariables2D *precalc, int i, int j)
{
    /*
    Calculates the right hand side of the elliptic equation in 2D.

    Parameters
    ----------
    bg : struct
        A pointer to the BackgroundVariables struct.
    fg : struct
        A pointer to the ForegroundVariables2D struct.
    grid_info : struct
        A pointer to the GridInfo2D struct.
    precalc : struct
        A pointer to the PrecalculatedVariables struct.
    i : int
        The index of the current cell in the z-direction.
    j : int
        The index of the current cell in the y-direction.
    
    Returns
    -------
    rhs : FLOAT_P
        The right hand side of the elliptic equation.
    */

    FLOAT_P rhs = 0.0; // This is the return value

    // Creating pointers to foreground arrays
    FLOAT_P **rho1 = fg->rho1;
    FLOAT_P **vy = fg->vy;
    FLOAT_P **vz = fg->vz;

    // Creating pointers to background arrays
    FLOAT_P *rho0 = bg->rho0;
    FLOAT_P *g = bg->g;
    FLOAT_P *grad_g = precalc->grad_g;
    FLOAT_P *grad_rho0 = precalc->grad_rho0;

    // Calculate the derivatives
    // First derivatives
    FLOAT_P dvy_dy = central_first_derivative_y_2D(vy, i, j, grid_info->ny, precalc->one_over_2dy);
    FLOAT_P dvy_dz = central_first_derivative_z_2D(vy, i, j, precalc->one_over_2dz);
    FLOAT_P dvz_dy = central_first_derivative_y_2D(vz, i, j, grid_info->ny, precalc->one_over_2dy);
    FLOAT_P dvz_dz = central_first_derivative_z_2D(vz, i, j, precalc->one_over_2dz);

    FLOAT_P drho1_dz = central_first_derivative_z_2D(rho1, i, j, precalc->one_over_2dz);

    // Second derivatives
    FLOAT_P dd_vy_ddy = central_second_derivative_y_2D(vy, i, j, grid_info->ny, precalc->one_over_dydy);
    FLOAT_P dd_vz_ddz = central_second_derivative_z_2D(vz, i, j, precalc->one_over_dzdz);

    // Mixed derivatives
    FLOAT_P dd_vy_dydz = central_second_derivative_yz_2D(vy, i, j, grid_info->ny, precalc->one_over_4dydz);
    FLOAT_P dd_vz_dydz = central_second_derivative_yz_2D(vz, i, j, grid_info->ny, precalc->one_over_4dydz);

    // Squared of derivatives
    FLOAT_P dvy_dy_sqrd = dvy_dy*dvy_dy;
    FLOAT_P dvz_dz_sqrd = dvz_dz*dvz_dz;

    #if GRAVITY_ON == 1
    {
        rhs -= g[i]*drho1_dz + rho1[i][j]*grad_g[i];
    }
    #endif // GRAVITY_ON

    #if ADVECTION_ON == 1
    {
        rhs -= rho0[i] * ( vy[i][j]*dd_vy_ddy + vz[i][j]*dd_vz_ddz + dvy_dy_sqrd + dvz_dz_sqrd
                          +2*dvz_dy*dvy_dz + vy[i][j]*dd_vz_dydz + vz[i][j]*dd_vy_dydz )
               +grad_rho0[i] * (vy[i][j]*dvz_dy + vz[i][j]*dvz_dz);
    }
    #endif // ADVECTION_ON

    #if VISCOSITY_ON == 1
        FLOAT_P ddd_vz_dydydz = central_third_derivative_yyz_2D(vz, i, j, grid_info->ny, precalc->one_over_8dydydz);
        FLOAT_P ddd_vy_dydzdz = central_third_derivative_yzz_2D(vy, i, j, grid_info->ny, precalc->one_over_8dydzdz);

        rhs += precalc->two_VIS_COEFF*(ddd_vy_dydzdz + ddd_vz_dydydz);
    #endif // VISCOSITY_ON
    
    return rhs;
}