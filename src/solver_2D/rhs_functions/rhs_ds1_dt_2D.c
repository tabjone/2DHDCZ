#include "rhs_functions.h"

FLOAT_P rhs_ds1_dt_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, int i, int j)
{
    /*
    Calculates the right hands side of the entropy equation in 2D.

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
    FLOAT_P *T0 = bg->T0;

    FLOAT_P ds1_dy = upwind_first_derivative_y(s1, vy, i, j, dy, ny);
    FLOAT_P ds1_dz = upwind_first_derivative_z(s1, vz, i, j, dz, nz_full);
    


    #if ADVECTION_ON == 1
        rhs -= vy[i][j]*ds1_dy + vz[i][j]*ds1_dz + vz[i][j]*grad_s0[i];
    #endif // ADVECTION_ON

    #if VISCOSITY_ON == 1
        FLOAT_P dvy_dz = central_first_derivative_z(vy, i, j, dz, nz_full);
        FLOAT_P dvz_dy = central_first_derivative_y(vz, i, j, dy, ny);
        
        rhs += VISCOSITY_COEFF/(bg->T0[i]*bg->rho0[i]) * (dvy_dz + dvz_dy)*(dvy_dz + dvz_dy);
    #endif // VISCOSITY_ON

    #if THERMAL_DIFFUSIVITY_ON == 1
        FLOAT_P dT0_dz = (T0[i+1] - T0[i-1])/(2*dz);
        FLOAT_P dd_s1_ddy = central_second_derivative_y(s1, i, j, dy, ny);
        FLOAT_P dd_s1_ddz = central_second_derivative_z(s1, i, j, dz, nz_full);

        rhs += 1.866e6/(bg->T0[i]*bg->rho0[i]) * (dT0_dz*(ds1_dy+ds1_dz) + bg->T0[i]*(dd_s1_ddy + dd_s1_ddz));
    #endif // THERMAL_DIFFUSIVITY_ON
    
    return rhs;
}