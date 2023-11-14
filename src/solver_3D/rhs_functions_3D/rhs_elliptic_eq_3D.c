#include "rhs_functions_3D.h"

FLOAT_P rhs_elliptic_eq_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, int i, int j, int k)
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
    /*
    // Getting the grid info
    int nx = grid_info->nx;
    int ny = grid_info->ny;
    int nz = grid_info->nz;
    FLOAT_P dx = grid_info->dx;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P dz = grid_info->dz;

    // Creating pointers to foreground arrays
    FLOAT_P ***rho1 = fg->rho1;
    FLOAT_P ***vx = fg->vx;
    FLOAT_P ***vy = fg->vy;
    FLOAT_P ***vz = fg->vz;

    // Creating pointers to background arrays
    FLOAT_P *rho0 = bg->rho0;
    FLOAT_P *g = bg->g;
    FLOAT_P *grad_g = bg->grad_g;
    FLOAT_P *grad_rho0 = bg->grad_rho0;

    // Periodic boundary conditions
    int j_minus = periodic_boundary(j-1, ny);
    int j_plus = periodic_boundary(j+1, ny);
    int k_minus = periodic_boundary(k-1, nx);
    int k_plus = periodic_boundary(k+1, nx);

    // Calculate the derivatives

    // Mixed derivatives
    FLOAT_P dd_vx_dxdz, dd_vz_dxdz, dd_vz_dydz, dd_vy_dydz, dd_vx_dxdy, dd_vy_dxdy;

    // Mixed derivatives
    dd_vy_dydz = (vy[i+1][j_plus][k]-vy[i+1][j_minus][k]-vy[i-1][j_plus][k]+vy[i-1][j_minus][k])/((4.0*dy*dz));
    dd_vz_dydz = (vz[i+1][j_plus][k]-vz[i+1][j_minus][k]-vz[i-1][j_plus][k]+vz[i-1][j_minus][k])/((4.0*dy*dz));
    dd_vz_dxdz = (vz[i+1][j][k_plus]-vz[i+1][j][k_minus]-vz[i-1][j][k_plus]+vz[i-1][j][k_minus])/((4.0*dx*dz));
    dd_vx_dxdz = (vx[i+1][j][k_plus]-vx[i+1][j][k_minus]-vx[i-1][j][k_plus]+vx[i-1][j][k_minus])/((4.0*dx*dz));
    dd_vx_dxdy = (vx[i][j_plus][k_plus]-vx[i][j_plus][k_minus]-vx[i][j_minus][k_plus]+vx[i][j_minus][k_minus])/((4.0*dx*dy));
    dd_vy_dxdy = (vy[i][j_plus][k_plus]-vy[i][j_plus][k_minus]-vy[i][j_minus][k_plus]+vy[i][j_minus][k_minus])/((4.0*dx*dy));


    // First derivatives
    FLOAT_P d_rho1_dz = central_first_derivative_z_3D(rho1, i, j, k, dz, nz);
    FLOAT_P dvx_dx = central_first_derivative_x_3D(vx, i, j, k, dx, nx);
    FLOAT_P dvx_dy = central_first_derivative_y_3D(vx, i, j, k, dy, ny);
    FLOAT_P dvx_dz = central_first_derivative_z_3D(vx, i, j, k, dz, nz);
    FLOAT_P dvy_dx = central_first_derivative_x_3D(vy, i, j, k, dx, nx);
    FLOAT_P dvy_dy = central_first_derivative_y_3D(vy, i, j, k, dy, ny);
    FLOAT_P dvy_dz = central_first_derivative_z_3D(vy, i, j, k, dz, nz);
    FLOAT_P dvz_dx = central_first_derivative_x_3D(vz, i, j, k, dx, nx);
    FLOAT_P dvz_dy = central_first_derivative_y_3D(vz, i, j, k, dy, ny);
    FLOAT_P dvz_dz = central_first_derivative_z_3D(vz, i, j, k, dz, nz);

    // Second derivatives
    FLOAT_P dd_vx_ddx = central_second_derivative_x_3D(vx, i, j, k, dx, nx);
    FLOAT_P dd_vy_ddy = central_second_derivative_y_3D(vy, i, j, k, dy, ny);
    FLOAT_P dd_vz_ddz = central_second_derivative_z_3D(vz, i, j, k, dz, nz);

    #if GRAVITY_ON == 1
    {
        rhs -= g[i]*d_rho1_dz + rho1[i][j][k]*grad_g[i];
    }
    #endif // GRAVITY_ON

    FLOAT_P dvx_dx_sqrd = dvx_dx*dvx_dx;
    FLOAT_P dvy_dy_sqrd = dvy_dy*dvy_dy;
    FLOAT_P dvz_dz_sqrd = dvz_dz*dvz_dz;

    #if ADVECTION_ON == 1
    {
        rhs -= rho0[i]*(vx[i][j][k]*dd_vx_ddx + vy[i][j][k]*dd_vy_ddy + vz[i][j][k]*dd_vz_ddz + dvx_dx_sqrd + dvy_dy_sqrd + dvz_dz_sqrd + 2*(dvy_dx*dvx_dy + dvz_dx*dvx_dz + dvz_dy*dvy_dz)
        + vx[i][j][k] * dd_vy_dxdy + vx[i][j][k] * dd_vz_dxdz + vy[i][j][k] * dd_vx_dxdy + vy[i][j][k] * dd_vz_dydz + vz[i][j][k] * dd_vx_dxdz + vz[i][j][k] * dd_vy_dydz)
        + grad_rho0[i]* (vx[i][j][k]*dvz_dx + vy[i][j][k]*dvz_dy + vz[i][j][k]*dvz_dz);
    }
    #endif // ADVECTION_ON

    #if VISCOSITY_ON == 1
        FLOAT_P ddd_vx_dxddy = ((vx[i][j_plus][k_plus] - 2.0 * vx[i][j][k_plus] + vx[i][j_minus][k_plus]) -
                                (vx[i][j_plus][k_minus] - 2.0 * vx[i][j][k_minus] + vx[i][j_minus][k_minus]) ) /
                               (2.0 * dx * dy * dy);

        FLOAT_P ddd_vx_dxddz = ((vx[i+1][j][k_plus] - 2.0 * vx[i][j][k_plus] + vx[i-1][j][k_plus]) -
                                (vx[i+1][j][k_minus] - 2.0 * vx[i][j][k_minus] + vx[i-1][j][k_minus])) /
                               (2.0 * dx * dz * dz);

        FLOAT_P ddd_vy_ddxdy = ((vy[i][j_plus][k_plus] - 2.0 * vy[i][j_plus][k] + vy[i][j_plus][k_minus]) -
                                (vy[i][j_minus][k_plus] - 2.0 * vy[i][j_minus][k] + vy[i][j_minus][k_minus])) /
                               (2.0 * dx * dx * dy);

        FLOAT_P ddd_vy_dyddz = ((vy[i+1][j_plus][k] - 2.0 * vy[i][j_plus][k] + vy[i-1][j_plus][k]) -
                                (vy[i+1][j_minus][k] - 2.0 * vy[i][j_minus][k] + vy[i-1][j_minus][k])) /
                               (2.0 * dy * dz * dz);

        FLOAT_P ddd_vz_ddxdz = ((vz[i+1][j][k_plus] - 2.0 * vz[i+1][j][k] + vz[i+1][j][k_minus]) -
                                (vz[i-1][j][k_plus] - 2.0 * vz[i-1][j][k] + vz[i-1][j][k_minus])) /
                               (2.0 * dx * dx * dz);

        FLOAT_P ddd_vz_ddydz = ((vz[i+1][j_plus][k] - 2.0 * vz[i+1][j][k] + vz[i+1][j_minus][k]) -
                                (vz[i-1][j_plus][k] - 2.0 * vz[i-1][j][k] + vz[i-1][j_minus][k])) /
                               (2.0 * dy * dy * dz);

        rhs += 2.0*VISCOSITY_COEFF*(ddd_vx_dxddy + ddd_vx_dxddz + ddd_vy_ddxdy + ddd_vy_dyddz + ddd_vz_ddxdz + ddd_vz_ddydz);
    #endif // VISCOSITY_ON
    */

    return rhs;
}