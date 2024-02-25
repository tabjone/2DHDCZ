#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "MPI_module/mpi_info_struct.h"

#include "global_parameters.h"
#include "global_constants.h"
#include "solver/boundary/boundary_2D/boundary_2D.h"
#include "array_utilities/array_memory_management/array_memory_management.h"

void equation_of_state_2D(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info)
{
    /*
    Calculates the foreground density from the equation of state.

    Parameters
    ----------
    fg : ForegroundVariables2D
        A pointer to the ForegroundVariables2D struct.
    bg : BackgroundVariables
        A pointer to the BackgroundVariables struct.
    grid_info : GridInfo2D
        A pointer to the GridInfo2D struct.
    */

    // Getting grid info
    int ny = grid_info->ny;
    int nz_full = grid_info->nz_full;

    FLOAT_P c_p = K_B / (MU * M_U) /(1.0-1.0/GAMMA);

    for (int i = 0; i < nz_full; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->rho1[i][j] = (-1.0/GAMMA * fg->p1[i][j]/bg->p0[i] - fg->s1[i][j]/c_p) * bg->rho0[i];
        }
    }
}