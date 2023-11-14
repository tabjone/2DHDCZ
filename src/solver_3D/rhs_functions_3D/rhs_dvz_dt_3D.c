#include "rhs_functions_3D.h"

FLOAT_P rhs_dvz_dt_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, int i, int j, int k)
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
    /*
    // Getting the grid info
    int nx = grid_info->nx;
    int ny = grid_info->ny;
    int nz = grid_info->nz;
    int nz_full = grid_info->nz_full;
    FLOAT_P dx = grid_info->dx;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P dz = grid_info->dz;
    
    // Creating pointers to foreground arrays
    FLOAT_P ***rho1 = fg->rho1;
    FLOAT_P ***p1 = fg->p1;
    FLOAT_P ***vx = fg->vx;
    FLOAT_P ***vy = fg->vy;
    FLOAT_P ***vz = fg->vz;

    // Creating pointers to background arrays
    FLOAT_P *one_over_rho0 = bg->one_over_rho0;
    FLOAT_P *g = bg->g;

    // Calculate the derivatives
    FLOAT_P dp1_dz = central_first_derivative_z_3D(p1, i, j, k, dz, nz);

    FLOAT_P dvz_dx = upwind_first_derivative_x_3D(vz, vx, i, j, k, dx, nx);
    FLOAT_P dvz_dy = upwind_first_derivative_y_3D(vz, vy, i, j, k, dy, ny);
    FLOAT_P dvz_dz = upwind_first_derivative_z_3D(vz, vz, i, j, k, dz, nz);

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
        // Periodic boundary conditions
        int j_minus = periodic_boundary(j-1, ny);
        int j_plus = periodic_boundary(j+1, ny);
        int k_minus = periodic_boundary(k-1, nx);
        int k_plus = periodic_boundary(k+1, nx);

        FLOAT_P dd_vz_ddx = central_second_derivative_x_3D(vz, i, j, k, dx, nx);
        FLOAT_P dd_vz_ddy = central_second_derivative_y_3D(vz, i, j, k, dy, ny);

        FLOAT_P dd_vx_dxdz = (vx[i+1][j][k_plus] - vx[i+1][j][k_minus] - vx[i-1][j][k_plus] + vx[i-1][j][k_minus])/(4.0*dx*dz);
        FLOAT_P dd_vy_dydz = (vy[i+1][j_plus][k] - vy[i+1][j_minus][k] - vy[i-1][j_plus][k] + vy[i-1][j_minus][k])/ (4.0*dy*dz);
        
        rhs += VISCOSITY_COEFF*one_over_rho0[i]*(dd_vx_dxdz + dd_vz_ddx + dd_vy_dydz + dd_vz_ddy);
    #endif // VISCOSITY_ON
    */
    return rhs;
}