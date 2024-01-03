#include "data_structures/grid_info/grid_info_1D/grid_info_struct_1D.h"
#include "data_structures/foreground_data/foreground_data_1D/foreground_variables_struct_1D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "MPI_module/MPI_module.h"
#include "solver/equation_of_state/equation_of_state_1D/equation_of_state_1D.h"
#include "solver/first_law_of_thermodynamics/first_law_of_thermodynamics_1D/first_law_of_thermodynamics_1D.h"
#include "global_initialization.h"
#include "global_boundary.h"
#include <math.h>
#include "global_constants.h"
#include "global_parameters.h"
#include "initialize_foreground_1D.h"

void initialize_foreground_sod_shock_vertical_1D(struct ForegroundVariables1D *fg, struct BackgroundVariables *bg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info)
{
    /*
    Initializes the foreground struct with a Sod Shock Tube test.

    Parameters
    ----------
    fg : ForegroundVariables1D
        A pointer to the ForegroundVariables1D struct.
    bg : BackgroundVariables
        A pointer to the BackgroundVariables struct.
    grid_info : GridInfo1D
        A pointer to the GridInfo1D struct.
    */

    // Getting grid info
    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;
    int nz = grid_info->nz;

    FLOAT_P dz = grid_info->dz;
    FLOAT_P z0 = grid_info->z0;

    initialize_foreground_zeros_1D(fg, grid_info);

    FLOAT_P z, tanh_value;
    FLOAT_P width = (R_END - R_START) * R_SUN * 0.05;
    FLOAT_P steepness = 1.0;
    FLOAT_P midpoint = (R_START + (R_END - R_START) / 2.0) * R_SUN;

    for (int i = nz_ghost; i < nz_full - nz_ghost; i++)
    {
        z = (i-nz_ghost) * dz + z0;
        
        tanh_value = tanh(steepness * (z - midpoint) / width);
        fg->p1[i] = 0.5 * ((UPPER_PRESSURE_BOUNDARY - LOWER_PRESSURE_BOUNDARY) * tanh_value + UPPER_PRESSURE_BOUNDARY + LOWER_PRESSURE_BOUNDARY);
        fg->s1[i] = 0.5 * ((UPPER_ENTROPY_BOUNDARY - LOWER_ENTROPY_BOUNDARY) * tanh_value + UPPER_ENTROPY_BOUNDARY + LOWER_ENTROPY_BOUNDARY);
    }
    communicate_1D_ghost_above_below(fg->p1, mpi_info, nz, nz_ghost);
    communicate_1D_ghost_above_below(fg->s1, mpi_info, nz, nz_ghost);

    first_law_of_thermodynamics_1D(fg, bg, grid_info);
    equation_of_state_1D(fg, bg, grid_info);
}