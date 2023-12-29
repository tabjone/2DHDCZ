#include "global_float_precision.h"
#include "global_parameters.h"
#include "spacial_derivatives_module/derivatives_3D/derivatives_3D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_3D/foreground_variables_struct_3D.h"
#include "data_structures/grid_info/grid_info_3D/grid_info_struct_3D.h"
#include "data_structures/precalculated_data/precalculated_data_3D/precalculated_data_3D.h"

FLOAT_P rhs_dvz_dt_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct PrecalculatedVariables3D *precalc, int i, int j, int k)
{
    /*
    Calculates the right hands side of the momentum equation for the z-direction in 3D.

    Parameters
    ----------
    bg : struct
        A pointer to the BackgroundVariables struct.
    fg : struct
        A pointer to the ForegroundVariables struct.
    grid_info : struct
        A pointer to the GridInfo struct.
    precalc : struct
        A pointer to the PrecalculatedVariables struct.
    i : int
        The index of the current cell in the z-direction.
    j : int
        The index of the current cell in the y-direction.
    k : int
        The index of the current cell in the x-direction.

    Returns
    -------
    rhs : FLOAT_P
        The right hand side of the momentum equation for the z-direction.
    */

    FLOAT_P rhs = 0.0; // This is the return value
    
    // Getting the grid info
    int nx = grid_info->nx;
    int ny = grid_info->ny;
    int nz = grid_info->nz;
    int nz_full = grid_info->nz_full;
    
    // Creating pointers to foreground arrays
    FLOAT_P ***rho1 = fg->rho1;
    FLOAT_P ***p1 = fg->p1;
    FLOAT_P ***vx = fg->vx;
    FLOAT_P ***vy = fg->vy;
    FLOAT_P ***vz = fg->vz;

    // Creating pointers to background arrays
    FLOAT_P *one_over_rho0 = precalc->one_over_rho0;
    FLOAT_P *g = bg->g;

    // Calculate the derivatives
    FLOAT_P dp1_dz = central_first_derivative_z_3D(p1, i, j, k, precalc->one_over_2dz);

    FLOAT_P dvz_dx = upwind_first_derivative_x_3D(vz, vx, i, j, k, nx, precalc->one_over_dx, precalc->one_over_2dx);
    FLOAT_P dvz_dy = upwind_first_derivative_y_3D(vz, vy, i, j, k, ny, precalc->one_over_dy, precalc->one_over_2dy);
    FLOAT_P dvz_dz = upwind_first_derivative_z_3D(vz, vz, i, j, k, precalc->one_over_dz, precalc->one_over_2dz);

    #if GAS_PRESSURE_ON == 1
        rhs -= one_over_rho0[i]* dp1_dz;
    #endif // GAS_PRESSURE_ON

    #if GRAVITY_ON == 1
        rhs -= one_over_rho0[i]*rho1[i][j][k] * g[i];
    #endif // GRAVITY_ON

    // Advective term
    #if ADVECTION_ON == 1
        rhs -= vx[i][j][k]*dvz_dx + vy[i][j][k]*dvz_dy + vz[i][j][k]*dvz_dz;
    #endif // ADVECTION_ON

    #if VISCOSITY_ON == 1
        FLOAT_P dd_vz_ddx = central_second_derivative_x_3D(vz, i, j, k, nx, precalc->one_over_dxdx);
        FLOAT_P dd_vz_ddy = central_second_derivative_y_3D(vz, i, j, k, ny, precalc->one_over_dydy);

        FLOAT_P dd_vx_dxdz = central_second_derivative_xz_3D(vx, i, j, k, nx, precalc->one_over_4dxdz);
        FLOAT_P dd_vy_dydz = central_second_derivative_yz_3D(vy, i, j, k, ny, precalc->one_over_4dydz);
        
        rhs += precalc->VIS_COEFF_over_rho0[i]*(dd_vx_dxdz + dd_vz_ddx + dd_vy_dydz + dd_vz_ddy);
    #endif // VISCOSITY_ON
    
    return rhs;
}