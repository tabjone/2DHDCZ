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

    // Getting the grid info
    int nx = grid_info->nx;
    int ny = grid_info->ny;
    int nz = grid_info->nz;
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

    return rhs;
}