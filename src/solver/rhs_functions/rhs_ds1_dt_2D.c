#include "rhs_functions.h"

#if DIMENSIONS == 2
FLOAT_P rhs_ds1_dt_2D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j)
{
    /*
    Calculates the right hands side of the entropy equation in 2D.

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
        The right hand side of the entropy equation.
    */

    FLOAT_P rhs = 0.0; // This is the return value

    // Getting the grid info
    int ny = grid_info->ny;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P dz = grid_info->dz;

    // Creating pointers to foreground arrays
    FLOAT_P **vy = fg->vy;
    FLOAT_P **vz = fg->vz;
    FLOAT_P **s1 = fg->s1;

    // Creating pointers to background arrays
    FLOAT_P *grad_s0 = bg->grad_s0;

    // Periodic boundary conditions
    int j_minus = periodic_boundary(j-1, ny);
    int j_plus = periodic_boundary(j+1, ny); 
    #if UPWIND_ORDER > 1
        int j_minus2 = periodic_boundary(j-2, ny);
        int j_plus2 = periodic_boundary(j+2, ny);
    #endif

    // Calculate the derivatives
    FLOAT_P ds1_dy, ds1_dz;

    #if UPWIND_ORDER == 1
        if (vy[i][j] >= 0)
        {
            ds1_dy = backward_first_derivative_first_order(s1[i][j], s1[i][j_minus], dy);
        }
        else
        {
            ds1_dy = forward_first_derivative_first_order(s1[i][j], s1[i][j_plus], dy);
        }
        if (vz[i][j] >= 0)
        {
            ds1_dz = backward_first_derivative_first_order(s1[i][j], s1[i-1][j], dz);
        }
        else
        {
            ds1_dz = forward_first_derivative_first_order(s1[i][j], s1[i+1][j], dz);
        }
    #elif UPWIND_ORDER == 2
        if (vy[i][j] >= 0)
        {
            ds1_dy = backward_first_derivative_second_order(s1[i][j], s1[i][j_minus], s1[i][j_minus2], dy);
        }
        else
        {
            ds1_dy = forward_first_derivative_second_order(s1[i][j], s1[i][j_plus], s1[i][j_plus2], dy);
        }
        if (vz[i][j] >= 0)
        {
            ds1_dz = backward_first_derivative_second_order(s1[i][j], s1[i-1][j], s1[i-2][j], dz);
        }
        else
        {
            ds1_dz = forward_first_derivative_second_order(s1[i][j], s1[i+1][j], s1[i+2][j], dz);
        }
    #endif // UPWIND_ORDER

    #if ADVECTION_ON == 1
        rhs -= vy[i][j]*ds1_dy + vz[i][j]*ds1_dz + vz[i][j]*grad_s0[i];
    #endif // ADVECTION_ON == 1

    #if BFIELD_ON == 1
        FLOAT_P *eta_over_four_pi_rho0_T0 = bg->eta_over_four_pi_rho0_T0;

        FLOAT_P **Bx = fg->Bx;
        FLOAT_P **By = fg->By;
        FLOAT_P **Bz = fg->Bz;

        FLOAT_P dBx_dy, dBy_dz, dBx_dz, dBz_dy;
        #if CENTRAL_ORDER == 2
            dBx_dy = central_first_derivative_second_order(Bx[i][j_minus], Bx[i][j_plus], dy);
            dBy_dz = central_first_derivative_second_order(By[i-1][j], By[i+1][j], dz);
            dBx_dz = central_first_derivative_second_order(Bx[i-1][j], Bx[i+1][j], dz);
            dBz_dy = central_first_derivative_second_order(Bz[i][j_minus], Bz[i][j_plus], dy);
        #endif // CENTRAL_ORDER == 2

        rhs -= eta_over_four_pi_rho0_T0[i] * ((dBz_dy - dBy_dz)*(dBz_dy - dBy_dz) + (dBx_dz)*(dBx_dz) + (dBx_dy)*(dBx_dy));
    #endif // BFIELD_ON == 1

    return rhs;
}
#endif // DIMENSIONS == 2