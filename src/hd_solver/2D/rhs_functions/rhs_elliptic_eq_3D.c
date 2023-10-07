#include "rhs_functions.h"

#include "rhs_functions.h"

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
    FLOAT_P dvx_dy, dvx_dz, dvy_dx, dvy_dz, dvz_dy, dvz_dx;
    FLOAT_P drho1_dz;

    // Mixed derivatives
    FLOAT_P dd_vx_dxdz, dd_vz_dxdz, dd_vz_dydz, dd_vy_dydz;

    // Second derivatives
    FLOAT_P dd_vz_ddz;

    #if CENTRAL_ORDER == 2
    // First derivatives
    dvx_dy = central_first_derivative_second_order(vx[i][j_minus][k], vx[i][j_plus][k], dy);
    dvx_dz = central_first_derivative_second_order(vx[i-1][j][k], vx[i+1][j][k], dz);
    dvy_dx = central_first_derivative_second_order(vy[i][j][k_minus], vy[i][j][k_plus], dx);
    dvy_dz = central_first_derivative_second_order(vy[i-1][j][k], vy[i+1][j][k], dz);
    dvz_dy = central_first_derivative_second_order(vz[i][j_minus][k], vz[i][j_plus][k], dy);
    dvz_dx = central_first_derivative_second_order(vz[i][j][k_minus], vz[i][j][k_plus], dx);
    drho1_dz = central_first_derivative_second_order(rho1[i-1][j], rho1[i+1][j], dz);

    // Mixed derivatives
    dd_vy_dydz = (vy[i+1][j_plus][k]-vy[i+1][j_minus][k]-vy[i-1][j_plus][k]+vy[i-1][j_minus][k])/((4.0*dy*dz));
    dd_vz_dydz = (vz[i+1][j_plus][k]-vz[i+1][j_minus][k]-vz[i-1][j_plus][k]+vz[i-1][j_minus][k])/((4.0*dy*dz));
    dd_vz_dxdz = (vz[i+1][j][k_plus]-vz[i+1][j][k_minus]-vz[i-1][j][k_plus]+vz[i-1][j][k_minus])/((4.0*dx*dz));
    dd_vx_dxdz = (vx[i+1][j][k_plus]-vx[i+1][j][k_minus]-vx[i-1][j][k_plus]+vx[i-1][j][k_minus])/((4.0*dx*dz));
    
    // Second derivatives
    dd_vz_ddz = central_second_derivative_second_order(vz[i][j], vz[i-1][j], vz[i+1][j], dz);
    #endif
    
    #if GRAVITY_ON == 1
    {
        rhs -= g[i]*drho1_dz + rho1[i][j][k]*grad_g[i];
    }
    #endif // GRAVITY_ON

    #if ADVECTION_ON == 1
    {
        rhs -= rho0[i]*(dvy_dx*dvx_dy + 2*dvz_dx*dvx_dz + dvx_dy*dvy_dx + 2*dvz_dy*dvy_dz + vz[i][j][k]*(dd_vx_dxdz+dd_vy_dy_dz) + vx[i][j][k]*dd_vz_dxdz + vy[i][j][k]*dd_vz_dydz + vz[i][j][k]*dd_vz_ddz) - grad_rho0[i]*(vx[i][j][k]*dvz_dx + vy[i][j][k]*dvz_dy);
    }
    #endif // ADVECTION_ON

    return rhs;
}