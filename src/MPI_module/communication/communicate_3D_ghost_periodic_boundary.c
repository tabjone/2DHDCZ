#include <stdbool.h>
#include <mpi.h>
#include "../mpi_info_struct.h"
#include "global_float_precision.h"

void communicate_3D_ghost_periodic_boundary(FLOAT_P ***array, struct MpiInfo *mpi_info, int nz, int nz_ghost, int ny, int nx) 
{
    /*
    Communicates the ghost cells of the first process to the last process and vice versa for a 3D array.

    Parameters
    ----------
    array : float pointer
        Pointer to the array to be communicated.
    mpi_info : struct MpiInfo pointer
        Pointer to the struct containing mpi info.
    nz : int
        Size of the actual cells in the array in the z-direction.
    nz_ghost : int
        Number of ghost cells.
    ny : int
        Size of the array in the y-direction.
    nx : int
        Size of the array in the x-direction.
    */

    // Package size
    int n_send = nz_ghost*ny*nx;

    // Getting mpi info
    int size = mpi_info->size;

    // If only one process, we need to copy the top rows to the bottom ghost cells and vice versa. And also make sure the boundary on the top is the same as the bottom.
    if (mpi_info->size == 1)
    {
        for (int i = 0; i < nz_ghost; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                for (int k = 0; k < nx; k++)
                {
                    array[i][j][k] = array[nz+i-1][j][k]; // Copying top rows to bottom ghost cells.
                    array[nz+nz_ghost+i][j][k] = array[nz_ghost+i+1][j][k]; // Copying bottom rows to top ghost cells.
                }
            }
        }
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                array[nz+nz_ghost][j][k] = array[nz_ghost][j][k]; // Making sure the boundary on the top is the same as the bottom.
            }
        }
        return; // No communication needed.
    }

    bool has_neighbor_above = mpi_info->has_neighbor_above;
    bool has_neighbor_below = mpi_info->has_neighbor_below;

    MPI_Status status;
    // Sending and recieving to/from ghost cells on the other side of the domain.
    if (!has_neighbor_above)
    {
        // Send to bottom
        MPI_Send(&array[nz-1][0][0], n_send, MPI_FLOAT_P, 0, 0, MPI_COMM_WORLD);
        // Recieve from bottom
        MPI_Recv(&array[nz+nz_ghost][0][0], n_send, MPI_FLOAT_P, 0, 1, MPI_COMM_WORLD, &status);
    }
    // Send and recieve to/from top
    if (!has_neighbor_below)
    {
        // Recieve from top
        MPI_Recv(&array[0][0][0], n_send, MPI_FLOAT_P, size-1, 0, MPI_COMM_WORLD, &status);
        // Send to top
        MPI_Send(&array[nz_ghost+1][0][0], n_send, MPI_FLOAT_P, size-1, 1, MPI_COMM_WORLD);
    }
    // Lastly we make sure the boundary on the top is the same as the bottom.
    if (!has_neighbor_below)
    {
        MPI_Send(&array[nz_ghost][0][0], nx*ny, MPI_FLOAT_P, size-1, 0, MPI_COMM_WORLD);
    }
    if (!has_neighbor_above)
    {
        MPI_Recv(&array[nz+nz_ghost-1][0][0], nx*ny, MPI_FLOAT_P, 0, 0, MPI_COMM_WORLD, &status);
    }
}