#include "one_time_step_3D.h"

FLOAT_P one_time_step_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg_prev, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt_last, bool first_timestep)
{   
    /*
    Calculates the foreground at the next timestep using Runge-Kutta methods.

    Parameters
    ----------
    bg : struct
        A pointer to the BackgroundVariables struct.
    fg_prev : struct
        A pointer to the ForegroundVariables struct at the previous timestep.
    fg : struct
        A pointer to the ForegroundVariables struct at the current timestep.
    grid_info : struct
        A pointer to the GridInfo struct.
    mpi_info : struct
        A pointer to the MpiInfo struct.
    dt_last : FLOAT_P
        The timestep used at the previous timestep.
    first_timestep : bool
        True if this is the first timestep, false otherwise.

    Returns
    -------
    dt : FLOAT_P
        The timestep used at the current timestep.
    */

    FLOAT_P dt;

    #if TIME_ORDER == 1
        dt = rk1_3D(bg, fg_prev, fg, grid_info, mpi_info, dt_last, first_timestep);

    #elif TIME_ORDER == 2
        dt = rk2_3D(bg, fg_prev, fg, grid_info, mpi_info, dt_last, first_timestep);

    #elif TIME_ORDER == 3
        dt = rk3_3D(bg, fg_prev, fg, grid_info, mpi_info, dt_last, first_timestep);
    #endif // TIME_ORDER

    return dt;
}