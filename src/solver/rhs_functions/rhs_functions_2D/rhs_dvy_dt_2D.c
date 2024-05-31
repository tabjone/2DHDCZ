#include "global_float_precision.h"
#include "global_parameters.h"
#include "global_constants.h"
#include "spacial_derivatives_module/derivatives_2D/derivatives_2D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "data_structures/precalculated_data/precalculated_data_2D/precalculated_data_struct_2D.h"

#include <mpi.h>
#include <math.h>

FLOAT_P rhs_dvy_dt_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct PrecalculatedVariables2D *precalc, int i, int j)
{
    /*
    Calculates the right hands side of the momentum equation for the y-direction in 2D.

    Parameters
    ----------
    bg : struct
        A pointer to the BackgroundVariables struct.
    fg : struct
        A pointer to the ForegroundVariables2D struct.
    grid_info : struct
        A pointer to the GridInfo2D struct.
    precalc : struct
        A pointer to the PrecalculatedVariables struct.
    i : int
        The index of the current cell in the z-direction.
    j : int
        The index of the current cell in the y-direction.

    Returns
    -------
    rhs : FLOAT_P
        The right hand side of the momentum equation for the y-direction.
    */

    FLOAT_P rhs = 0.0; // This is the return value
    
    // Creating pointers to foreground arrays
    FLOAT_P **p1 = fg->p1;
    FLOAT_P **vy = fg->vy;
    FLOAT_P **vz = fg->vz;

    // Creating pointers to background arrays
    FLOAT_P *one_over_rho0 = precalc->one_over_rho0;

    // Calculate the derivatives
    FLOAT_P dp1_dy = central_first_derivative_y_2D(p1, i, j, grid_info->ny, precalc->one_over_2dy);
    FLOAT_P dvy_dy = upwind_first_derivative_y_2D(vy, vy, i, j, grid_info->ny, precalc->one_over_dy, precalc->one_over_2dy);
    FLOAT_P dvy_dz = upwind_first_derivative_z_2D(vy, vz, i, j, precalc->one_over_dz, precalc->one_over_2dz);
 
    #if GAS_PRESSURE_ON == 1
        #if COORDINATES == 0
            rhs -= one_over_rho0[i]* dp1_dy;
        #elif COORDINATES == 1
            rhs -= one_over_rho0[i]/bg->r[i] * dp1_dy;
        #endif // COORDINATES
    #endif // GAS_PRESSURE_ON

    #if CORIOLIS_ON == 1
        FLOAT_P omega, f;
        int n = 2;

        f = (OMEGA_CORE - OMEGA_EQ) * pow(bg->r[i]/R_SUN, n);

        omega = OMEGA_EQ + f;
        rhs += 2.0*omega*vz[i][j];
    #endif // CORIOLIS_ON

    // Advective term
    #if ADVECTION_ON == 1
        #if COORDINATES == 0
            rhs -= vy[i][j]*dvy_dy + vz[i][j]*dvy_dz;
        #elif COORDINATES == 1
            rhs -= vz[i][j]*dvy_dz + vy[i][j]/bg->r[i]*dvy_dy;
        #endif // COORDINATES
    #endif // ADVECTION_ON

    #if VISCOSITY_ON == 1
        FLOAT_P dd_vy_ddy = central_second_derivative_y_2D(vy, i, j, grid_info->ny, precalc->one_over_dydy);
        FLOAT_P dd_vz_dydz = central_second_derivative_yz_2D(vz, i, j, grid_info->ny, precalc->one_over_4dydz);
        FLOAT_P dd_vy_ddz = central_second_derivative_z_2D(vy, i, j, precalc->one_over_dzdz);        

        rhs += precalc->VIS_COEFF_over_rho0[i]*
            (
                4.0/3.0*dd_vy_ddy + 1.0/3.0 * dd_vz_dydz + dd_vy_ddz
            );
    #endif // VISCOSITY_ON

    return rhs;
}