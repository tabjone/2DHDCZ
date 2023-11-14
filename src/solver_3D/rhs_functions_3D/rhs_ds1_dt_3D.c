#include "rhs_functions_3D.h"

FLOAT_P rhs_ds1_dt_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, int i, int j, int k)
{
    /*
    Calculates the right hands side of the entropy equation in 3D.

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
        The right hand side of the entropy equation.
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
    FLOAT_P ***vx = fg->vx;
    FLOAT_P ***vy = fg->vy;
    FLOAT_P ***vz = fg->vz;
    FLOAT_P ***s1 = fg->s1;

    // Creating pointers to background arrays
    FLOAT_P *grad_s0 = bg->grad_s0;

    // Calculate the derivatives
    FLOAT_P ds1_dx, ds1_dy, ds1_dz;

    ds1_dx = upwind_first_derivative_x_3D(s1, vx, i, j, k, dx, nx);
    ds1_dy = upwind_first_derivative_y_3D(s1, vy, i, j, k, dy, ny);
    ds1_dz = upwind_first_derivative_z_3D(s1, vz, i, j, k, dz, nz);

    #if ADVECTION_ON == 1
        rhs -= vx[i][j][k]*ds1_dx + vy[i][j][k]*ds1_dy + vz[i][j][k]*(ds1_dz + grad_s0[i]);
    #endif // ADVECTION_ON
    
    #if THERMAL_DIFFUSIVITY_ON == 1
        FLOAT_P *T0 = bg->T0;

        FLOAT_P dT0_dz = (T0[i+1] - T0[i-1])/(2*dz);
        FLOAT_P dd_s1_ddx = central_second_derivative_x_3D(s1, i, j, k, dx, nx);
        FLOAT_P dd_s1_ddy = central_second_derivative_y_3D(s1, i, j, k, dy, ny);
        FLOAT_P dd_s1_ddz = central_second_derivative_z_3D(s1, i, j, k, dz, nz_full);

        rhs += 1.866e6/(bg->T0[i]*bg->rho0[i]) * (dT0_dz*(ds1_dx + ds1_dy+ds1_dz) + bg->T0[i]*(dd_s1_ddx + dd_s1_ddy + dd_s1_ddz));
    #endif // THERMAL_DIFFUSIVITY_ON

    #if VISCOSITY_ON == 1
        FLOAT_P dvy_dz = central_first_derivative_z_3D(vy, i, j, k, dz, nz_full);
        FLOAT_P dvz_dy = central_first_derivative_y_3D(vz, i, j, k, dy, ny);
        FLOAT_P dvx_dz = central_first_derivative_z_3D(vx, i, j, k, dz, nz_full);
        FLOAT_P dvz_dx = central_first_derivative_x_3D(vz, i, j, k, dx, nx);
        FLOAT_P dvx_dy = central_first_derivative_y_3D(vx, i, j, k, dy, ny);
        FLOAT_P dvy_dx = central_first_derivative_x_3D(vy, i, j, k, dx, nx);
        
        rhs += VISCOSITY_COEFF/(bg->T0[i]*bg->rho0[i]) * ((dvy_dz + dvz_dy)*(dvy_dz + dvz_dy) + (dvx_dy + dvy_dx)*(dvx_dy + dvy_dx) + (dvx_dz + dvz_dx)*(dvx_dz + dvz_dx));
    #endif // VISCOSITY_ON
    */
    return rhs;
}