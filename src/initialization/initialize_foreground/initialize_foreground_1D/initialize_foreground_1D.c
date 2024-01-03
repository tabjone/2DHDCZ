#include "initialize_foreground_1D.h"
#include "global_initialization.h"

void initialize_foreground_1D(struct ForegroundVariables1D *fg, struct BackgroundVariables *bg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info)
{
    /*
    Initializes the foreground struct.

    Parameters
    ----------
    fg : ForegroundVariables1D
        A pointer to the ForegroundVariables1D struct.
    bg : BackgroundVariables
        A pointer to the BackgroundVariables struct.
    grid_info : GridInfo1D
        A pointer to the GridInfo1D struct.
    mpi_info : MpiInfo
        A pointer to the MpiInfo struct.
    */

    #if IC_ENTROPY_PERTURBATION == 1
        initialize_foreground_entropy_perturbations_1D(fg, bg, grid_info, mpi_info);
    #elif IC_SOD_SHOCK == 1
        initialize_foreground_sod_shock_vertical_1D(fg, bg, grid_info, mpi_info);
    #elif IC_OSCILLATION_MODES == 1
        initialize_foreground_oscillation_modes_1D(fg, bg, grid_info, mpi_info);
    #elif IC_ZEROS == 1
        initialize_foreground_zeros_1D(fg, grid_info);
    #elif IC_RANDOM_OSCILLATIONS == 1
        initialize_foreground_random_oscillations_1D(fg, bg, grid_info, mpi_info);
    #endif // IC_TYPE

}
