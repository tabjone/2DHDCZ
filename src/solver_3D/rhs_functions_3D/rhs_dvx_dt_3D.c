#include "rhs_functions_3D.h"

FLOAT_P rhs_dvx_dt_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, int i, int j, int k)
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
    FLOAT_P ***p1 = fg->p1;
    FLOAT_P ***vx = fg->vx;
    FLOAT_P ***vy = fg->vy;
    FLOAT_P ***vz = fg->vz;

    // Creating pointers to background arrays
    FLOAT_P *one_over_rho0 = bg->one_over_rho0;

    // Calculate the derivatives
    FLOAT_P dp1_dx = central_first_derivative_x_3D(p1, i, j, k, dx, nx);

    FLOAT_P dvx_dx = upwind_first_derivative_x_3D(vx, vx, i, j, k, dx, nx);
    FLOAT_P dvx_dy = upwind_first_derivative_y_3D(vx, vy, i, j, k, dy, ny);
    FLOAT_P dvx_dz = upwind_first_derivative_z_3D(vx, vz, i, j, k, dz, nz);

    #if GAS_PRESSURE_ON == 1
        rhs -= one_over_rho0[i]* dp1_dx;
    #endif // GAS_PRESSURE_ON

    // Advective term
    #if ADVECTION_ON == 1
        rhs -= vx[i][j][k]*dvx_dx + vy[i][j][k]*dvx_dy + vz[i][j][k]*dvx_dz;
    #endif // ADVECTION_ON

    #if VISCOSITY_ON == 1
        // Periodic boundary conditions
        int j_minus = periodic_boundary(j-1, ny);
        int j_plus = periodic_boundary(j+1, ny);
        int k_minus = periodic_boundary(k-1, nx);
        int k_plus = periodic_boundary(k+1, nx);

        FLOAT_P dd_vx_ddy = central_second_derivative_y_3D(vx, i, j, k, dy, ny);
        FLOAT_P dd_vx_ddz = central_second_derivative_z_3D(vx, i, j, k, dz, nz_full);

        FLOAT_P dd_vy_dxdy = (vy[i][j_plus][k_plus] - vy[i][j_minus][k_minus] - vy[i][j_minus][k_plus] + vy[i][j_minus][k_minus])/(4.0*dx*dy);
        FLOAT_P dd_vz_dxdz = (vz[i+1][j][k_plus] - vz[i+1][j][k_minus] - vz[i-1][j][k_plus] + vz[i-1][j][k_minus])/(4.0*dx*dz);

        rhs += VISCOSITY_COEFF*one_over_rho0[i]*(dd_vx_ddy + dd_vy_dxdy + dd_vx_ddz + dd_vz_dxdz);
    #endif // VISCOSITY_ON
    */
    return rhs;
}