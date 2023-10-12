#include <stdlib.h>
#include "initialization.h"
#include "global_parameters.h"

void initialize_foreground_random(struct ForegroundVariables *fg, struct BackgroundVariables *bg , struct GridInfo *grid_info)
{
    /*
    Initializes the foreground struct with random values.

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

    initialize_foreground_zeros(fg, grid_info); // Start by setting everyhing to zero

    // Getting grid info
    int ny = grid_info->ny;
    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;

    // Inside the grid we pick random values for p1 and s1 and let the boundaries stay zero
    for (int i = nz_ghost; i < nz_full-nz_ghost-1; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            fg->p1[i][j] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 1e-5;  // random between -1e-5 to 1e-5
            fg->s1[i][j] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 20000;  // random between -20000 to 20000
        }
    }

    // Solving first law of thermodynamics and equation of state to find rho1 and T1
    first_law_thermodynamics(fg, bg, grid_info);
    equation_of_state(fg, bg, grid_info);
}
