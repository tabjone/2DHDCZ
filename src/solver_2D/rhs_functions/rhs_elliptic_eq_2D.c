#include "rhs_functions.h"

FLOAT_P rhs_elliptic_eq_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, int i, int j)
{
    /*
    Calculates the right hand side of the elliptic equation in 2D.

    Parameters
    ----------
    bg : struct
        A pointer to the BackgroundVariables struct.
    fg : struct
        A pointer to the ForegroundVariables2D struct.
    grid_info : struct
        A pointer to the GridInfo2D struct.
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
    int nz_full = grid_info->nz_full;
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

    // Calculate the derivatives
    // First derivatives
    FLOAT_P dvy_dy = central_first_derivative_y(vy, i, j, dy, ny);
    FLOAT_P dvy_dz = central_first_derivative_z(vy, i, j, dz, nz_full);
    FLOAT_P dvz_dy = central_first_derivative_y(vz, i, j, dy, ny);
    FLOAT_P dvz_dz = central_first_derivative_z(vz, i, j, dz, nz_full);

    FLOAT_P drho1_dz = central_first_derivative_z(rho1, i, j, dz, nz_full);

    // Second derivatives
    FLOAT_P dd_vy_ddy = central_second_derivative_y(vy, i, j, dy, ny);
    FLOAT_P dd_vz_ddz = central_second_derivative_z(vz, i, j, dz, nz_full);

    // Mixed derivatives
    FLOAT_P dd_vy_dydz = (vy[i+1][j_plus] - vy[i+1][j_minus] - vy[i-1][j_plus] + vy[i-1][j_minus])/(4.0*dy*dz);
    FLOAT_P dd_vz_dydz = (vz[i+1][j_plus] - vz[i+1][j_minus] - vz[i-1][j_plus] + vz[i-1][j_minus])/(4.0*dy*dz);
    
    // Squared of derivatives
    FLOAT_P dvy_dy_sqrd = dvy_dy*dvy_dy;
    FLOAT_P dvz_dz_sqrd = dvz_dz*dvz_dz;

    #if GRAVITY_ON == 1
    {
        rhs -= g[i]*drho1_dz + rho1[i][j]*grad_g[i];
    }
    #endif // GRAVITY_ON

    #if ADVECTION_ON == 1
    {
        rhs -= rho0[i]*(vy[i][j]*dd_vy_ddy + vz[i][j]*dd_vz_ddz + dvy_dy_sqrd + dvz_dz_sqrd + 2*dvz_dy*dvy_dz + vy[i][j]*dd_vz_dydz + vz[i][j]*dd_vy_dydz)
        + grad_rho0[i] * (vy[i][j]*dvz_dy + vz[i][j]*dvz_dz);
    }
    #endif // ADVECTION_ON

    #if VISCOSITY_ON == 1
        FLOAT_P dd_vz_ddy = central_second_derivative_y(vz, i, j, dy, ny);

        FLOAT_P ddd_vy_dyddz = (vy[i+1][j_plus] - vy[i+1][j_minus] - 2.0*vy[i][j_plus] + 2.0*vy[i][j_minus] + vy[i-1][j_plus] - vy[i-1][j_minus])/(2.0*dy*dz*dz);


        rhs += 2*VISCOSITY_COEFF*(ddd_vy_dyddz + dd_vz_ddy);
    #endif // VISCOSITY_ON
    
    return rhs;
}