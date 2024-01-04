#include "spacial_derivatives_module/derivatives_3D/derivatives_3D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_3D/foreground_variables_struct_3D.h"
#include "data_structures/grid_info/grid_info_3D/grid_info_struct_3D.h"
#include "data_structures/precalculated_data/precalculated_data_3D/precalculated_data_3D.h"
#include "global_float_precision.h"
#include "global_parameters.h"

FLOAT_P rhs_dvx_dt_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct PrecalculatedVariables3D *precalc, int i, int j, int k)
{
    /*
    Calculates the right hands side of the momentum equation for the x-direction in 3D.

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
        The right hand side of the momentum equation for the x-direction.
    */

    FLOAT_P rhs = 0.0; // This is the return value
    
    // Getting the grid info
    int nx = grid_info->nx;
    int ny = grid_info->ny;
    
    // Creating pointers to foreground arrays
    FLOAT_P ***p1 = fg->p1;
    FLOAT_P ***vx = fg->vx;
    FLOAT_P ***vy = fg->vy;
    FLOAT_P ***vz = fg->vz;

    // Creating pointers to background arrays
    FLOAT_P *one_over_rho0 = precalc->one_over_rho0;

    // Calculate the derivatives
    FLOAT_P dp1_dx = central_first_derivative_x_3D(p1, i, j, k, nx, precalc->one_over_2dx);

    FLOAT_P dvx_dx = upwind_first_derivative_x_3D(vx, vx, i, j, k, nx, precalc->one_over_dx, precalc->one_over_2dx);
    FLOAT_P dvx_dy = upwind_first_derivative_y_3D(vx, vy, i, j, k, ny, precalc->one_over_dy, precalc->one_over_2dy);
    FLOAT_P dvx_dz = upwind_first_derivative_z_3D(vx, vz, i, j, k, precalc->one_over_dz, precalc->one_over_2dz);

    #if GAS_PRESSURE_ON == 1
        rhs -= one_over_rho0[i]* dp1_dx;
    #endif // GAS_PRESSURE_ON

    // Advective term
    #if ADVECTION_ON == 1
        rhs -= vx[i][j][k]*dvx_dx + vy[i][j][k]*dvx_dy + vz[i][j][k]*dvx_dz;
    #endif // ADVECTION_ON

    #if VISCOSITY_ON == 1
        FLOAT_P dd_vx_ddy = central_second_derivative_y_3D(vx, i, j, k, ny, precalc->one_over_dydy);
        FLOAT_P dd_vx_ddz = central_second_derivative_z_3D(vx, i, j, k, precalc->one_over_dzdz);

        FLOAT_P dd_vy_dxdy = central_second_derivative_xy_3D(vy, i, j, k, nx, ny, precalc->one_over_4dxdy);
        FLOAT_P dd_vz_dxdz = central_second_derivative_xz_3D(vz, i, j, k, nx, precalc->one_over_4dxdz);

        rhs += precalc->VIS_COEFF_over_rho0[i]*(dd_vx_ddy + dd_vy_dxdy + dd_vx_ddz + dd_vz_dxdz);
    #endif // VISCOSITY_ON
    
    return rhs;
}