#include "rhs_functions.h"
#include "global_parameters.h"

#if DIMENSIONS == 3
FLOAT_P rhs_elliptic_eq_3D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j, int k)
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

    // Getting the grid info
    int nx = grid_info->nx;
    int ny = grid_info->ny;
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
    // First derivatives
    FLOAT_P d_vx_dy, d_vx_dz, d_vy_dx, d_vy_dz, d_vz_dy, d_vz_dx, d_vx_dx, d_vy_dy, d_vz_dz;
    FLOAT_P d_rho1_dz;

    // Mixed derivatives
    FLOAT_P dd_vx_dxdz, dd_vz_dxdz, dd_vz_dydz, dd_vy_dydz, dd_vx_dxdy, dd_vy_dxdy;

    // Second derivatives
    FLOAT_P dd_vz_ddz, dd_vx_ddx, dd_vy_ddy;

    #if CENTRAL_ORDER == 2
        // First derivatives
        d_vx_dy = central_first_derivative_second_order(vx[i][j_minus][k], vx[i][j_plus][k], dy);
        d_vx_dz = central_first_derivative_second_order(vx[i-1][j][k], vx[i+1][j][k], dz);
        d_vy_dx = central_first_derivative_second_order(vy[i][j][k_minus], vy[i][j][k_plus], dx);
        d_vy_dz = central_first_derivative_second_order(vy[i-1][j][k], vy[i+1][j][k], dz);
        d_vz_dy = central_first_derivative_second_order(vz[i][j_minus][k], vz[i][j_plus][k], dy);
        d_vz_dx = central_first_derivative_second_order(vz[i][j][k_minus], vz[i][j][k_plus], dx);
        d_vx_dx = central_first_derivative_second_order(vx[i][j][k_minus], vx[i][j][k_plus], dx);
        d_vy_dy = central_first_derivative_second_order(vy[i][j_minus][k], vy[i][j_plus][k], dy);
        d_vz_dz = central_first_derivative_second_order(vz[i-1][j][k], vz[i+1][j][k], dz);
        d_rho1_dz = central_first_derivative_second_order(rho1[i-1][j][k], rho1[i+1][j][k], dz);

        // Mixed derivatives
        dd_vy_dydz = (vy[i+1][j_plus][k]-vy[i+1][j_minus][k]-vy[i-1][j_plus][k]+vy[i-1][j_minus][k])/((4.0*dy*dz));
        dd_vz_dydz = (vz[i+1][j_plus][k]-vz[i+1][j_minus][k]-vz[i-1][j_plus][k]+vz[i-1][j_minus][k])/((4.0*dy*dz));
        dd_vz_dxdz = (vz[i+1][j][k_plus]-vz[i+1][j][k_minus]-vz[i-1][j][k_plus]+vz[i-1][j][k_minus])/((4.0*dx*dz));
        dd_vx_dxdz = (vx[i+1][j][k_plus]-vx[i+1][j][k_minus]-vx[i-1][j][k_plus]+vx[i-1][j][k_minus])/((4.0*dx*dz));
        dd_vx_dxdy = (vx[i][j_plus][k_plus]-vx[i][j_plus][k_minus]-vx[i][j_minus][k_plus]+vx[i][j_minus][k_minus])/((4.0*dx*dy));
        dd_vy_dxdy = (vy[i][j_plus][k_plus]-vy[i][j_plus][k_minus]-vy[i][j_minus][k_plus]+vy[i][j_minus][k_minus])/((4.0*dx*dy));
        
        // Second derivatives
        dd_vz_ddz = central_second_derivative_second_order(vz[i][j][k], vz[i-1][j][k], vz[i+1][j][k], dz);
        dd_vx_ddx = central_second_derivative_second_order(vx[i][j][k], vx[i][j][k_minus], vx[i][j][k_plus], dx);
        dd_vy_ddy = central_second_derivative_second_order(vy[i][j][k], vy[i][j_minus][k], vy[i][j_plus][k], dy);
    #endif
    
    #if GRAVITY_ON == 1
    {
        rhs -= g[i]*drho1_dz + rho1[i][j][k]*grad_g[i];
    }
    #endif // GRAVITY_ON

    #if ADVECTION_ON == 1
    {
        rhs -= rho0[i]*(vx[i][j][k]*dd_vx_ddx + vy[i][j][k]*dd_vy_ddy + vz[i][j][k]*dd_vz_ddz + (d_vx_dx)*(d_vx_dx) + (d_vy_dy)*(d_vy_dy) + (d_vz_dz)*(d_vz_dz) + 2*(d_vy_dx*d_vx_dy + d_vz_dx*d_vx_dz + d_vz_dy*d_vy_dz)
        + vx[i][j][k] * dd_vy_dxdy + vx[i][j][k] * dd_vz_dxdz + vy[i][j][k] * dd_vx_dxdy + vy[i][j][k] * dd_vz_dydz + vz[i][j][k] * dd_vx_dxdz + vz[i][j][k] * dd_vy_dydz)
        - d_rho1_dz* (vx[i][j][k]*d_vx_dz + vy[i][j][k]*d_vy_dz + vz[i][j][k]*d_vz_dz);
    }
    #endif // ADVECTION_ON

    return rhs;
}
#endif // DIMENSIONS == 3