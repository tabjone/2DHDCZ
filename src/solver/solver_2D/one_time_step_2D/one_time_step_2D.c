#include "../../spacial_derivatives_module/spacial_derivatives_module.h"
#include "global_float_precision.h"
#include "../../data_structures/background_data/background_variables_struct.h"
#include "../../data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "../../data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "one_time_step.h"

FLOAT_P one_time_step(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables *precalc, FLOAT_P dt_last, bool first_timestep)
{   
    /*
    Calculates the foreground at the next timestep using Runge-Kutta methods.

    Parameters
    ----------
    bg : struct
        A pointer to the BackgroundVariables struct.
    fg_prev : struct
        A pointer to the ForegroundVariables2D struct at the previous timestep.
    fg : struct
        A pointer to the ForegroundVariables2D struct at the current timestep.
    grid_info : struct
        A pointer to the GridInfo2D struct.
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
        dt = rk1_2D(bg, fg_prev, fg, grid_info, mpi_info, precalc, dt_last, first_timestep);
    #elif TIME_ORDER == 2
        dt = rk2_2D(bg, fg_prev, fg, grid_info, mpi_info, precalc, dt_last, first_timestep);
    #elif TIME_ORDER == 3
        dt = rk3_2D(bg, fg_prev, fg, grid_info, mpi_info, precalc, dt_last, first_timestep);
    #endif // TIME_ORDER

    return dt;
}