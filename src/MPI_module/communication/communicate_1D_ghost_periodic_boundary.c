#include <stdbool.h>
#include <mpi.h>
#include "../mpi_info_struct.h"
#include "global_float_precision.h"

void communicate_1D_ghost_periodic_boundary(FLOAT_P *array, struct MpiInfo *mpi_info, int nz, int nz_ghost)
{
    /*
    Communicates the ghost cells of the first process to the last process and vice versa for a 1D array.

    Parameters
    ----------
    array : float pointer
        Pointer to the array to be communicated.
    mpi_info : struct MpiInfo pointer
        Pointer to the struct containing mpi info.
    nz : int
        Size of the actual cells in the array.
    nz_ghost : int
        Number of ghost cells.
    */
    
    // Getting mpi info
    int size = mpi_info->size;

    // If only one process, we need to copy the top rows to the bottom ghost cells and vice versa. And also make sure the boundary on the top is the same as the bottom.
    if (mpi_info->size == 1)
    {
        for (int i = 0; i < nz_ghost; i++)
        {
            array[i] = array[nz+i-1]; // Copying top rows to bottom ghost cells.
            array[nz+nz_ghost+i] = array[nz_ghost+i+1]; // Copying bottom rows to top ghost cells.
        }
            array[nz+nz_ghost] = array[nz_ghost]; // Making sure the boundary on the top is the same as the bottom.
        return; // No communication needed.
    }

    bool has_neighbor_above = mpi_info->has_neighbor_above;
    bool has_neighbor_below = mpi_info->has_neighbor_below;

    MPI_Status status;
    // Sending and recieving to/from ghost cells on the other side of the domain.
    if (!has_neighbor_above)
    {
        // Send to bottom
        MPI_Send(&array[nz-1], nz_ghost, MPI_FLOAT_P, 0, 0, MPI_COMM_WORLD);
        // Recieve from bottom
        MPI_Recv(&array[nz+nz_ghost], nz_ghost, MPI_FLOAT_P, 0, 1, MPI_COMM_WORLD, &status);
    }
    // Send and recieve to/from top
    if (!has_neighbor_below)
    {
        // Recieve from top
        MPI_Recv(&array[0], nz_ghost, MPI_FLOAT_P, size-1, 0, MPI_COMM_WORLD, &status);
        // Send to top
        MPI_Send(&array[nz_ghost+1], nz_ghost, MPI_FLOAT_P, size-1, 1, MPI_COMM_WORLD);
    }
    // Lastly we make sure the boundary on the top is the same as the bottom.
    if (!has_neighbor_below)
    {
        MPI_Send(&array[nz_ghost], 1, MPI_FLOAT_P, size-1, 0, MPI_COMM_WORLD);
    }
    if (!has_neighbor_above)
    {
        MPI_Recv(&array[nz+nz_ghost-1], 1, MPI_FLOAT_P, 0, 0, MPI_COMM_WORLD, &status);
    }
}