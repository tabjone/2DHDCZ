#include "initialization_2D.h"

void initialize_foreground(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info)
{
    /*
    Initializes the foreground struct.

    Parameters
    ----------
    fg : ForegroundVariables2D
        A pointer to the ForegroundVariables2D struct.
    bg : BackgroundVariables
        A pointer to the BackgroundVariables struct.
    grid_info : GridInfo2D
        A pointer to the GridInfo2D struct.
    mpi_info : MpiInfo
        A pointer to the MpiInfo struct.
    */
    #if MPI_ON == 0
        #if IC_ENTROPY_PERTUBATION == 1
            initialize_foreground_entropy_pertubation(fg, bg, grid_info, mpi_info);
        #elif IC_DENSITY_PERTUBATION == 1
            initialize_foreground_density_pertubation(fg, bg, grid_info);
        #elif IC_SOD_SHOCK == 1
            #if IC_SOD_SHOCK_DIRECTION == 0
                 sod_shock_horizontal(fg, bg, grid_info);
            #elif IC_SOD_SHOCK_DIRECTION == 1
                sod_shock_vertical(fg, bg, grid_info);
            #endif // IC_SOD_SHOCK_DIRECTION
        #endif // IC_TYPE

    #elif MPI_ON == 1
        #if IC_ENTROPY_PERTUBATION == 1
            initialize_foreground_entropy_pertubation_mpi(fg, bg, grid_info, mpi_info);
        #elif IC_DENSITY_PERTUBATION == 1
            //initialize_foreground_density_pertubation_mpi(fg, bg, grid_info);
        #elif IC_SOD_SHOCK == 1
            #if IC_SOD_SHOCK_DIRECTION == 0
                 //sod_shock_horizontal_mpi(fg, bg, grid_info);
            #elif IC_SOD_SHOCK_DIRECTION == 1
                //sod_shock_vertical_mpi(fg, bg, grid_info);
            #endif // IC_SOD_SHOCK_DIRECTION
        #endif // IC_TYPE
    #endif // MPI_ON
}
