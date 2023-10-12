#include "rhs_functions.h"

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
    int nz_full = grid_info->nz_full;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P dz = grid_info->dz;

    // Creating pointers to foreground arrays
    FLOAT_P **vy = fg->vy;
    FLOAT_P **vz = fg->vz;
    FLOAT_P **s1 = fg->s1;

    // Creating pointers to background arrays
    FLOAT_P *grad_s0 = bg->grad_s0;

    FLOAT_P ds1_dy = upwind_first_derivative_y(s1, vy, i, j, dy, ny);
    FLOAT_P ds1_dz = upwind_first_derivative_z(s1, vz, i, j, dz, nz_full);

    #if ADVECTION_ON == 1
        rhs -= vy[i][j]*ds1_dy + vz[i][j]*ds1_dz + vz[i][j]*grad_s0[i];
    #endif // ADVECTION_ON == 1

    return rhs;
}