#include "rhs_functions.h"

#if DIMENSIONS == 2
FLOAT_P rhs_elliptic_eq_2D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j)
{
    /*
    Calculates the right hand side of the elliptic equation in 2D.

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
    
    Returns
    -------
    rhs : FLOAT_P
        The right hand side of the elliptic equation.
    */

    FLOAT_P rhs = 0.0; // This is the return value

    // Getting the grid info
    int ny = grid_info->ny;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P dz = grid_info->dz;

    // Creating pointers to foreground arrays
    FLOAT_P **rho1 = fg->rho1;
    FLOAT_P **vy = fg->vy;
    FLOAT_P **vz = fg->vz;

    // Creating pointers to background arrays
    FLOAT_P *rho0 = bg->rho0;
    FLOAT_P *g = bg->g;
    FLOAT_P *grad_g = bg->grad_g;
    FLOAT_P *grad_rho0 = bg->grad_rho0;


    // Periodic boundary conditions
    int j_minus = periodic_boundary(j-1, ny);
    int j_plus = periodic_boundary(j+1, ny);

    // Write out terms for the mixed derivatives, u for up, d for down, r for right, l for left
    FLOAT_P vy_ur = vy[i+1][j_plus];
    FLOAT_P vy_ul = vy[i+1][j_minus];
    FLOAT_P vy_dr = vy[i-1][j_plus];
    FLOAT_P vy_dl = vy[i-1][j_minus];

    FLOAT_P vz_ur = vz[i+1][j_plus];
    FLOAT_P vz_ul = vz[i+1][j_minus];
    FLOAT_P vz_dr = vz[i-1][j_plus];
    FLOAT_P vz_dl = vz[i-1][j_minus];

    // Calculate the derivatives
    // First derivatives
    FLOAT_P d_vy_dy, d_vy_dz, d_vz_dy, d_vz_dz;
    FLOAT_P d_rho1_dz;
    // Second derivatives
    FLOAT_P dd_vy_ddy, dd_vz_ddz;
    // Mixed derivatives
    FLOAT_P dd_vy_dydz, dd_vz_dydz;

    #if CENTRAL_ORDER == 2
        // First derivatives
        d_vy_dy = central_first_derivative_second_order(vy[i][j_minus], vy[i][j_plus], dy);
        d_vy_dz = central_first_derivative_second_order(vy[i-1][j], vy[i+1][j], dz);
        d_vz_dy = central_first_derivative_second_order(vz[i][j_minus], vz[i][j_plus], dy);
        d_vz_dz = central_first_derivative_second_order(vz[i-1][j], vz[i+1][j], dz);
        d_rho1_dz = central_first_derivative_second_order(rho1[i-1][j], rho1[i+1][j], dz);

        // Second derivatives
        dd_vy_ddy = central_second_derivative_second_order(vy[i][j], vy[i][j_minus], vy[i][j_plus], dy);
        dd_vz_ddz = central_second_derivative_second_order(vz[i][j], vz[i-1][j], vz[i+1][j], dz);

        // Mixed derivatives
        dd_vy_dydz = (vy_ur-vy_ul-vy_dr+vy_dl)/(4.0*dy*dz);
        dd_vz_dydz = (vz_ur-vz_ul-vz_dr+vz_dl)/(4.0*dy*dz);
    #endif

    #if GRAVITY_ON == 1
    {
        rhs -= g[i]*drho1_dz + rho1[i][j]*grad_g[i];
    }
    #endif // GRAVITY_ON

    #if ADVECTION_ON == 1
    {
        rhs -= rho0[i]*(vy[i][j]*dd_vy_ddy + vz[i][j]*dd_vz_ddz + (d_vy_dy)*(d_vy_dy) + (d_vz_dz)*(d_vz_dz) + 2*d_vz_dy*d_vy_dz + vy[i][j]*dd_vz_dydz + vz[i][j]*dd_vy_dydz)
        + d_rho1_dz * (vy[i][j]*d_vz_dy + vz[i][j]*d_vz_dz);
    }
    #endif // ADVECTION_ON

    return rhs;
}
#endif // DIMENSIONS == 2