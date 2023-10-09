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
        bool has_neighbour_above = mpi_info->has_neighbour_above;
        bool has_neighbour_below = mpi_info->has_neighbour_below;
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
            // We first create random numbers for the boundary if there is a neighbour above or below
            if (has_neighbour_above)
            {
                fg->p1[nz_full-nz_ghost-1] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 1e-5;  // random between -1e-5 to 1e-5
                fg->s1[nz_full-nz_ghost-1] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 20000;  // random between -20000 to 20000
            }
            if (has_neighbour_below)
            {
                fg->p1[nz_ghost] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 1e-5;  // random between -1e-5 to 1e-5
                fg->s1[nz_ghost] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 20000;  // random between -20000 to 20000
            }

            // Getting mpi info
            bool has_neighbour_above = mpi_info->has_neighbour_above;
            bool has_neighbour_below = mpi_info->has_neighbour_below;
            int rank = mpi_info->rank;
            int size = mpi_info->size;

            // Then we communicate the boundary values to the neighbouring ranks
            MPI_Status status;
            // If rank is even, send first, then receive
            if (rank % 2 == 0) 
            {
                if (has_neighbour_above) 
                {
                    // Send to the next rank (above)
                    MPI_Send(&(fg->p1[nz_full-2*nz_ghost]), nz_ghost, MPI_FLOAT_P, rank + 1, 0, MPI_COMM_WORLD);
                    MPI_Send(&(fg->s1[nz_full-2*nz_ghost]), nz_ghost, MPI_FLOAT_P, rank + 1, 1, MPI_COMM_WORLD);
                }
                if (has_neighbour_below) 
                {
                    // Receive from the previous rank (below)
                    MPI_Recv(&(fg->p1[0]), nz_ghost, MPI_FLOAT_P, rank - 1, 0, MPI_COMM_WORLD, &status);
                    MPI_Recv(&(fg->s1[0]), nz_ghost, MPI_FLOAT_P, rank - 1, 1, MPI_COMM_WORLD, &status);
                }
            }
            // If rank is odd, receive first, then send
            else 
            {
                if (has_neighbour_above) {
                    // Receive from the previous rank (below)
                    MPI_Recv(&(fg->p1[0]), nz_ghost, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &status);
                    MPI_Recv(&(fg->s1[0]), nz_ghost, MPI_FLOAT, rank - 1, 1, MPI_COMM_WORLD, &status);
                }
                if (has_neighbour_below) 
                {
                    // Send to the next rank (above)
                    MPI_Send(&(fg->p1[nz_full-2*nz_ghost]), nz_ghost, MPI_FLOAT_P, rank + 1, 0, MPI_COMM_WORLD);
                    MPI_Send(&(fg->s1[nz_full-2*nz_ghost]), nz_ghost, MPI_FLOAT_P, rank + 1, 1, MPI_COMM_WORLD);
                }
            }
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
            // We first create random numbers for the boundary if there is a neighbour above or below
            if (has_neighbour_above)
            {
                for (int j = 0; j < ny; j++)
                {
                    fg->p1[nz_full-nz_ghost-1][j] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 1e-5;  // random between -1e-5 to 1e-5
                    fg->s1[nz_full-nz_ghost-1][j] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 20000;  // random between -20000 to 20000
                }
            }
            if (has_neighbour_below)
            {
                for (int j = 0; j < ny; j++)
                {
                    fg->p1[nz_ghost][j] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 1e-5;  // random between -1e-5 to 1e-5
                    fg->s1[nz_ghost][j] = (2.0 * (rand() / (FLOAT_P)RAND_MAX) - 1.0) * 20000;  // random between -20000 to 20000
                }
            }

            // Then we communicate the boundary values to the neighbouring ranks
            MPI_Status status;
            // If rank is even, send first, then receive
            if (rank % 2 == 0) 
            {
                if (has_neighbour_above) 
                {
                    // Send to the next rank (above)
                    MPI_Send(&(fg->p1[nz_full-2*nz_ghost][0]), nz_ghost*ny, MPI_FLOAT_P, rank + 1, 0, MPI_COMM_WORLD);
                    MPI_Send(&(fg->s1[nz_full-2*nz_ghost][0]), nz_ghost*ny, MPI_FLOAT_P, rank + 1, 1, MPI_COMM_WORLD);
                }
                if (has_neighbour_below) 
                {
                    // Receive from the previous rank (below)
                    MPI_Recv(&(fg->p1[0][0]), nz_ghost*ny, MPI_FLOAT_P, rank - 1, 0, MPI_COMM_WORLD, &status);
                    MPI_Recv(&(fg->s1[0][0]), nz_ghost*ny, MPI_FLOAT_P, rank - 1, 1, MPI_COMM_WORLD, &status);
                }
            }
            // If rank is odd, receive first, then send
            else 
            {
                if (has_neighbour_above) 
                {
                    // Receive from the previous rank (below)
                    MPI_Recv(&(fg->p1[0][0]), nz_ghost*ny, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &status);
                    MPI_Recv(&(fg->s1[0][0]), nz_ghost*ny, MPI_FLOAT, rank - 1, 1, MPI_COMM_WORLD, &status);
                }
                if (has_neighbour_below) 
                {
                    // Send to the next rank (above)
                    MPI_Send(&(fg->p1[nz_full-2*nz_ghost][0]), nz_ghost*ny, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD);
                    MPI_Send(&(fg->s1[nz_full-2*nz_ghost][0]), nz_ghost*ny, MPI_FLOAT, rank + 1, 1, MPI_COMM_WORLD);
                }
            }
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
            // We first create random numbers for the boundary if there is a neighbour above or below
            if (has_neighbour_above)
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
            if (has_neighbour_below)
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

            // Getting mpi info
            bool has_neighbour_above = mpi_info->has_neighbour_above;
            bool has_neighbour_below = mpi_info->has_neighbour_below;
            int rank = mpi_info->rank;
            int size = mpi_info->size;

            // Then we communicate the boundary values to the neighbouring ranks
            MPI_Status status;
            // If rank is even, send first, then receive
            if (rank % 2 == 0) 
            {
                if (has_neighbour_above) 
                {
                    // Send to the next rank (above)
                    MPI_Send(&(fg->p1[nz_full-2*nz_ghost][0][0]), nz_ghost*ny*nx, MPI_FLOAT_P, rank + 1, 0, MPI_COMM_WORLD);
                    MPI_Send(&(fg->s1[nz_full-2*nz_ghost][0][0]), nz_ghost*ny*nx, MPI_FLOAT_P, rank + 1, 1, MPI_COMM_WORLD);
                }
                if (has_neighbour_below) 
                {
                    // Receive from the previous rank (below)
                    MPI_Recv(&(fg->p1[0][0][0]), nz_ghost*ny*nx, MPI_FLOAT_P, rank - 1, 0, MPI_COMM_WORLD, &status);
                    MPI_Recv(&(fg->s1[0][0][0]), nz_ghost*ny*nx, MPI_FLOAT_P, rank - 1, 1, MPI_COMM_WORLD, &status);
                }
            }
            // If rank is odd, receive first, then send
            else 
            {
                if (has_neighbour_above) 
                {
                    // Receive from the previous rank (below)
                    MPI_Recv(&(fg->p1[0][0][0]), nz_ghost*ny*nx, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &status);
                    MPI_Recv(&(fg->s1[0][0][0]), nz_ghost*ny*nx, MPI_FLOAT, rank - 1, 1, MPI_COMM_WORLD, &status);
                }
                if (has_neighbour_below) 
                {
                    // Send to the next rank (above)
                    MPI_Send(&(fg->p1[nz_full-2*nz_ghost][0][0]), nz_ghost*ny*nx, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD);
                    MPI_Send(&(fg->s1[nz_full-2*nz_ghost][0][0]), nz_ghost*ny*nx, MPI_FLOAT, rank + 1, 1, MPI_COMM_WORLD);
                }
            }
        #endif // MPI_ON

        // Solving first law of thermodynamics and equation of state to find rho1 and T1
        first_law_thermodynamics(fg, bg, grid_info);
        equation_of_state(fg, bg, grid_info);
    #endif // DIMENSIONS
}
