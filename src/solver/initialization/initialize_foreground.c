#include "initialization.h"

void initialize_foreground(struct ForegroundVariables *fg, struct BackgroundVariables *bg, struct GridInfo *grid_info, struct MpiInfo *mpi_info)
{
    /*
    Initializes the foreground struct.

    Parameters
    ----------
    fg : ForegroundVariables
        A pointer to the ForegroundVariables struct.
    bg : BackgroundVariables
        A pointer to the BackgroundVariables struct.
    grid_info : GridInfo
        A pointer to the GridInfo struct.
    mpi_info : MpiInfo
        A pointer to the MpiInfo struct.
    */

    #if INITIALIZATION_TYPE == 0
        initialize_foreground_zeros(fg, grid_info);
    #elif INITIALIZATION_TYPE == 1
        initialize_foreground_velocity_right(fg, grid_info, mpi_info);
    #elif INITIALIZATION_TYPE == 2
        initialize_foreground_density_pertubation(fg, bg, grid_info, mpi_info);
    #elif INITIALIZATION_TYPE == 3
        initialize_foreground_random(fg, bg, grid_info);
    #endif // INITIALIZATION_TYPE
}
