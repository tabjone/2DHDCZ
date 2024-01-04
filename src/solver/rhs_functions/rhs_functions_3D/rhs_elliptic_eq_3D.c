#include "global_float_precision.h"
#include "global_parameters.h"
#include "spacial_derivatives_module/derivatives_3D/derivatives_3D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_3D/foreground_variables_struct_3D.h"
#include "data_structures/grid_info/grid_info_3D/grid_info_struct_3D.h"
#include "data_structures/precalculated_data/precalculated_data_3D/precalculated_data_3D.h"

FLOAT_P rhs_elliptic_eq_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct PrecalculatedVariables3D *precalc, int i, int j, int k)
{
    /*
    Calculates the right hand side of the elliptic equation in 3D.

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
        The right hand side of the elliptic equation.
    */

    FLOAT_P rhs = 0.0; // This is the return value
    
    // Getting the grid info
    int nx = grid_info->nx;
    int ny = grid_info->ny;

    // Creating pointers to foreground arrays
    FLOAT_P ***rho1 = fg->rho1;
    FLOAT_P ***vx = fg->vx;
    FLOAT_P ***vy = fg->vy;
    FLOAT_P ***vz = fg->vz;

    // Creating pointers to background arrays
    FLOAT_P *rho0 = bg->rho0;
    FLOAT_P *g = bg->g;
    FLOAT_P *grad_g = precalc->grad_g;
    FLOAT_P *grad_rho0 = precalc->grad_rho0;

    // Mixed derivatives
    FLOAT_P dd_vx_dxdz, dd_vz_dxdz, dd_vz_dydz, dd_vy_dydz, dd_vx_dxdy, dd_vy_dxdy;

    dd_vy_dydz = central_second_derivative_yz_3D(vy, i, j, k, ny, precalc->one_over_4dydz);
    dd_vz_dydz = central_second_derivative_yz_3D(vz, i, j, k, ny, precalc->one_over_4dydz);
    dd_vz_dxdz = central_second_derivative_xz_3D(vz, i, j, k, nx, precalc->one_over_4dxdz);
    dd_vx_dxdz = central_second_derivative_xz_3D(vx, i, j, k, nx, precalc->one_over_4dxdz);
    dd_vx_dxdy = central_second_derivative_xy_3D(vx, i, j, k, nx, ny, precalc->one_over_4dxdy);
    dd_vy_dxdy = central_second_derivative_xy_3D(vy, i, j, k, nx, ny, precalc->one_over_4dxdy);

    // First derivatives
    FLOAT_P d_rho1_dz = central_first_derivative_z_3D(rho1, i, j, k, precalc->one_over_2dz);
    FLOAT_P dvx_dx = central_first_derivative_x_3D(vx, i, j, k, nx, precalc->one_over_2dx);
    FLOAT_P dvx_dy = central_first_derivative_y_3D(vx, i, j, k, ny, precalc->one_over_2dy);
    FLOAT_P dvx_dz = central_first_derivative_z_3D(vx, i, j, k, precalc->one_over_2dz);
    FLOAT_P dvy_dx = central_first_derivative_x_3D(vy, i, j, k, nx, precalc->one_over_2dx);
    FLOAT_P dvy_dy = central_first_derivative_y_3D(vy, i, j, k, ny, precalc->one_over_2dy);
    FLOAT_P dvy_dz = central_first_derivative_z_3D(vy, i, j, k, precalc->one_over_2dz);
    FLOAT_P dvz_dx = central_first_derivative_x_3D(vz, i, j, k, nx, precalc->one_over_2dx);
    FLOAT_P dvz_dy = central_first_derivative_y_3D(vz, i, j, k, ny, precalc->one_over_2dy);
    FLOAT_P dvz_dz = central_first_derivative_z_3D(vz, i, j, k, precalc->one_over_2dz);

    // Second derivatives
    FLOAT_P dd_vx_ddx = central_second_derivative_x_3D(vx, i, j, k, nx, precalc->one_over_dxdx);
    FLOAT_P dd_vy_ddy = central_second_derivative_y_3D(vy, i, j, k, ny, precalc->one_over_dydy);
    FLOAT_P dd_vz_ddz = central_second_derivative_z_3D(vz, i, j, k, precalc->one_over_dzdz);

    #if GRAVITY_ON == 1
        rhs -= g[i]*d_rho1_dz + rho1[i][j][k]*grad_g[i];
    #endif // GRAVITY_ON

    FLOAT_P dvx_dx_sqrd = dvx_dx*dvx_dx;
    FLOAT_P dvy_dy_sqrd = dvy_dy*dvy_dy;
    FLOAT_P dvz_dz_sqrd = dvz_dz*dvz_dz;

    #if ADVECTION_ON == 1
        rhs -= rho0[i]*(vx[i][j][k]*dd_vx_ddx + vy[i][j][k]*dd_vy_ddy + vz[i][j][k]*dd_vz_ddz + dvx_dx_sqrd + dvy_dy_sqrd + dvz_dz_sqrd + 2*(dvy_dx*dvx_dy + dvz_dx*dvx_dz + dvz_dy*dvy_dz)
        + vx[i][j][k] * dd_vy_dxdy + vx[i][j][k] * dd_vz_dxdz + vy[i][j][k] * dd_vx_dxdy + vy[i][j][k] * dd_vz_dydz + vz[i][j][k] * dd_vx_dxdz + vz[i][j][k] * dd_vy_dydz)
        + grad_rho0[i]* (vx[i][j][k]*dvz_dx + vy[i][j][k]*dvz_dy + vz[i][j][k]*dvz_dz);
    #endif // ADVECTION_ON

    #if VISCOSITY_ON == 1
        FLOAT_P ddd_vx_dxdydy = central_third_derivative_xyy_3D(vx, i, j, k, nx, ny, precalc->one_over_8dxdydy);
        FLOAT_P ddd_vx_dxdzdz = central_third_derivative_xzz_3D(vx, i, j, k, nx, precalc->one_over_8dxdzdz);
        FLOAT_P ddd_vy_dxdxdy = central_third_derivative_xxy_3D(vy, i, j, k, nx, ny, precalc->one_over_8dxdxdy);
        FLOAT_P ddd_vy_dydzdz = central_third_derivative_yzz_3D(vy, i, j, k, ny, precalc->one_over_8dydzdz);
        FLOAT_P ddd_vz_dxdxdz = central_third_derivative_xxz_3D(vz, i, j, k, nx, precalc->one_over_8dxdxdz);
        FLOAT_P ddd_vz_dydydz = central_third_derivative_yyz_3D(vz, i, j, k, ny, precalc->one_over_8dydydz);

        rhs += precalc->two_VIS_COEFF*(ddd_vx_dxdydy + ddd_vx_dxdzdz + ddd_vy_dxdxdy + ddd_vy_dydzdz + ddd_vz_dxdxdz + ddd_vz_dydydz);
    #endif // VISCOSITY_ON
    
    return rhs;
}