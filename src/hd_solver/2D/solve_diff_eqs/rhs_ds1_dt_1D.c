#include "rhs_functions.h"

FLOAT_P rhs_ds1_dt_1D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i)
{
    /*
    Calculates the right hands side of the entropy equation in 1D.

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
    
    Returns
    -------
    rhs : FLOAT_P
        The right hand side of the entropy equation.
    */

    FLOAT_P rhs = 0.0; // This is the return value

    // Getting the grid info
    FLOAT_P dz = grid_info->dz;

    // Creating pointers to foreground arrays
    FLOAT_P *vy = fg->vy;
    FLOAT_P *vz = fg->vz;
    FLOAT_P *s1 = fg->s1;

    // Creating pointers to background arrays
    FLOAT_P *grad_s0 = bg->grad_s0;

    // Calculate the derivatives
    FLOAT_P ds1_dz;
    #if UPWIND_ORDER == 1
        if (vz[i] >= 0)
        {
            ds1_dz = backward_first_derivative_first_order(s1[i], s1[i-1], dz);
        }
        else
        {
            ds1_dz = forward_first_derivative_first_order(s1[i], s1[i+1], dz);
        }
    #elif UPWIND_ORDER == 2
        if (vz[i] >= 0)
        {
            ds1_dz = backward_first_derivative_second_order(s1[i], s1[i-1], s1[i-2], dz);
        }
        else
        {
            ds1_dz = forward_first_derivative_second_order(s1[i], s1[i+1], s1[i+2], dz);
        }
    #endif

    #if ADVECTION_ON == 1
        rhs -= vz[i]*ds1_dz;
    #endif

    return rhs;
}