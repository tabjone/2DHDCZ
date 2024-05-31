#include "initialize_foreground_2D.h"
#include "global_initialization.h"
#include "global_parameters.h"
#include "global_constants.h"
#include "MPI_module/MPI_module.h"
#include "solver/equation_of_state/equation_of_state_2D/equation_of_state_2D.h"

void initialize_foreground_2D(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info)
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

    #if IC_ENTROPY_PERTURBATION == 1
        initialize_foreground_entropy_perturbations_2D(fg, bg, grid_info, mpi_info);
    #elif IC_SOD_SHOCK == 1
        #if IC_SOD_SHOCK_DIRECTION == 0
            initialize_foreground_sod_shock_horizontal_2D(fg, bg, grid_info, mpi_info);
        #elif IC_SOD_SHOCK_DIRECTION == 1
            initialize_foreground_sod_shock_vertical_2D(fg, bg, grid_info, mpi_info);
        #endif // IC_SOD_SHOCK_DIRECTION
    #elif IC_OSCILLATION_MODES == 1
        initialize_foreground_oscillation_modes_2D(fg, bg, grid_info, mpi_info);
    #elif IC_ZEROS == 1
        initialize_foreground_zeros_2D(fg, grid_info);
    #elif IC_RANDOM_OSCILLATIONS == 1
        initialize_foreground_random_oscillations_2D(fg, bg, grid_info, mpi_info);
    #endif // IC_TYPE
    /*
    int nz_full = grid_info->nz_full;
    int nz = grid_info->nz;
    int nz_ghost = grid_info->nz_ghost;
    int ny = grid_info->ny;
    // velocity to the right with waves
    for (int i = 0; i < nz_full; i++)
        for (int j = 0; j < ny; j++)
            fg->vy[i][j] = 1.0e5;

    FLOAT_P Ly = Y_SIZE * R_SUN;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P y;
    FLOAT_P frequency = 2.0;
    FLOAT_P c_p = K_B / (MU * M_U) / (1.0 - 1.0/GAMMA);
    for (int i = 0; i < nz_full; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            y = j * dy;
            fg->s1[i][j] = sin(2.0 * M_PI * y / Ly * frequency) * 0.01 * c_p;
            fg->p1[i][j] = -bg->p0[i] * fg->s1[i][j]/c_p;
        }
    }

    communicate_2D_ghost_above_below(fg->s1, mpi_info, nz, nz_ghost, ny);
    communicate_2D_ghost_above_below(fg->p1, mpi_info, nz, nz_ghost, ny);

    equation_of_state_2D(fg, bg, grid_info, mpi_info); 
    */
}
