#include <stdlib.h>
#include "initialization.h"
#include "global_parameters.h"

void initialize_foreground_random(struct ForegroundVariables *fg, struct BackgroundVariables *bg , struct GridInfo *grid_info, struct MpiInfo *mpi_info)
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

    #if MPI_ON == 1
        // Getting mpi info
        bool has_neighbor_above = mpi_info->has_neighbor_above;
        bool has_neighbor_below = mpi_info->has_neighbor_below;
        int rank = mpi_info->rank;
        int size = mpi_info->size;
    #endif // MPI_ON

    initialize_foreground_zeros(fg, grid_info); // Start by setting everyhing to zero

    #if DIMENSIONS == 1
        // Getting grid info
        int nz_full = grid_info->nz_full;
        int nz_ghost = grid_info->nz_ghost;

        // Inside the grid we pick random values for p1 and s1 and let the boundaries stay zero
        for (int i = nz_ghost; i < nz_full-nz_ghost-1; i++)
        {
            fg->p1[i] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 1e-5;  // random between -1e-5 to 1e-5
            fg->s1[i] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 20000;  // random between -20000 to 20000
        }

        #if MPI_ON == 1
            // We first create random numbers for the boundary if there is a neighbor above or below
            if (has_neighbor_above)
            {
                fg->p1[nz_full-nz_ghost-1] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 1e-5;  // random between -1e-5 to 1e-5
                fg->s1[nz_full-nz_ghost-1] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 20000;  // random between -20000 to 20000
            }
            if (has_neighbor_below)
            {
                fg->p1[nz_ghost] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 1e-5;  // random between -1e-5 to 1e-5
                fg->s1[nz_ghost] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 20000;  // random between -20000 to 20000
            }
            // Communicating ghost cells
            communicate_1D_ghost_above_below(fg->p1, grid_info, mpi_info);
            communicate_1D_ghost_above_below(fg->s1, grid_info, mpi_info);
        #endif // MPI_ON

        // Solving first law of thermodynamics and equation of state to find rho1 and T1
        first_law_thermodynamics(fg, bg, grid_info);
        equation_of_state(fg, bg, grid_info);
    
    #elif DIMENSIONS == 2
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

        #if MPI_ON == 1
            // We first create random numbers for the boundary if there is a neighbor above or below
            if (has_neighbor_above)
            {
                for (int j = 0; j < ny; j++)
                {
                    fg->p1[nz_full-nz_ghost-1][j] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 1e-5;  // random between -1e-5 to 1e-5
                    fg->s1[nz_full-nz_ghost-1][j] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 20000;  // random between -20000 to 20000
                }
            }
            if (has_neighbor_below)
            {
                for (int j = 0; j < ny; j++)
                {
                    fg->p1[nz_ghost][j] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 1e-5;  // random between -1e-5 to 1e-5
                    fg->s1[nz_ghost][j] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 20000;  // random between -20000 to 20000
                }
            }
            
            // Communicating ghost cells
            communicate_2D_ghost_above_below(fg->p1, grid_info, mpi_info);
            communicate_2D_ghost_above_below(fg->s1, grid_info, mpi_info);
        #endif // MPI_ON

        // Solving first law of thermodynamics and equation of state to find rho1 and T1
        first_law_thermodynamics(fg, bg, grid_info);
        equation_of_state(fg, bg, grid_info);

    #elif DIMENSIONS == 3
        // Inside the grid we pick random values for p1 and s1 and let the boundaries stay zero
        for (int i = nz_ghost; i < nz_full-nz_ghost-1; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                for (int k = 0; k < nx; k++)
                {
                    fg->p1[i][j][k] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 1e-5;  // random between -1e-5 to 1e-5
                    fg->s1[i][j][k] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 20000;  // random between -20000 to 20000
                }
            }
        }

        #if MPI_ON == 1
            // We first create random numbers for the boundary if there is a neighbor above or below
            if (has_neighbor_above)
            {
                for (int j = 0; j < ny; j++)
                {
                    for (int k = 0; k < nx; k++)
                    {
                        fg->p1[nz_full-nz_ghost-1][j][k] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 1e-5;  // random between -1e-5 to 1e-5
                        fg->s1[nz_full-nz_ghost-1][j][k] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 20000;  // random between -20000 to 20000
                    }
                }
            }
            if (has_neighbor_below)
            {
                for (int j = 0; j < ny; j++)
                {
                    for (int k = 0; k < nx; k++)
                    {
                        fg->p1[nz_ghost][j][k] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 1e-5;  // random between -1e-5 to 1e-5
                        fg->s1[nz_ghost][j][k] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 20000;  // random between -20000 to 20000
                    }
                }
            }

            // Communicating ghost cells
            communicate_3D_ghost_above_below(fg->p1, grid_info, mpi_info);
            communicate_3D_ghost_above_below(fg->s1, grid_info, mpi_info);
        #endif // MPI_ON

        // Solving first law of thermodynamics and equation of state to find rho1 and T1
        first_law_thermodynamics(fg, bg, grid_info);
        equation_of_state(fg, bg, grid_info);
    #endif // DIMENSIONS
}
