#include "initialize_foreground_3D.h"
#include "global_initialization.h"

void initialize_foreground_3D(struct ForegroundVariables3D *fg, struct BackgroundVariables *bg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info)
{
    /*
    Initializes the foreground struct.

    Parameters
    ----------
    fg : ForegroundVariables3D
        A pointer to the ForegroundVariables3D struct.
    bg : BackgroundVariables
        A pointer to the BackgroundVariables struct.
    grid_info : GridInfo3D
        A pointer to the GridInfo3D struct.
    mpi_info : MpiInfo
        A pointer to the MpiInfo struct.
    */

    #if IC_ENTROPY_PERTURBATION == 1
        initialize_foreground_entropy_perturbations_3D(fg, bg, grid_info, mpi_info);
    #elif IC_SOD_SHOCK == 1
        #if IC_SOD_SHOCK_DIRECTION == 0
            initialize_foreground_sod_shock_horizontal_3D(fg, bg, grid_info, mpi_info);
        #elif IC_SOD_SHOCK_DIRECTION == 1
            initialize_foreground_sod_shock_vertical_3D(fg, bg, grid_info, mpi_info);
        #endif // IC_SOD_SHOCK_DIRECTION
    #elif IC_OSCILLATION_MODES == 1
        initialize_foreground_oscillation_modes_3D(fg, bg, grid_info, mpi_info);
    #elif IC_ZEROS == 1
        initialize_foreground_zeros_3D(fg, grid_info);
    #elif IC_RANDOM_OSCILLATIONS == 1
        initialize_foreground_random_oscillations_3D(fg, bg, grid_info, mpi_info);
    #endif // IC_TYPE

}
