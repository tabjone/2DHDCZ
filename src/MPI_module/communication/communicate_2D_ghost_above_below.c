#include <mpi.h>
#include "global_boundary.h"
#include "communication.h"

void communicate_2D_ghost_above_below(FLOAT_P **array, struct MpiInfo *mpi_info, int nz, int nz_ghost, int ny) 
{
    /*
    Communicates ghost cells above and below the current process for a 2D array using the odd-even method.

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
    */

    // Getting mpi info
    int rank = mpi_info->rank;
    bool has_neighbor_above = mpi_info->has_neighbor_above;
    bool has_neighbor_below = mpi_info->has_neighbor_below;

    // Returning if only one process no periodic boundary conditions.
    #if PERIODIC_BOUNDARY_VERTICAL == 0
        if (mpi_info->size == 1) 
        {
            return; // No communication needed.
        }
    #endif

    // Package size
    int n_send = nz_ghost*ny;

    // Handle the in-between processes.
    MPI_Status status;
    // Sending top rows to bottom ghost cells.
    if (rank % 2 == 0) {
        if (has_neighbor_above) {
            MPI_Send(&array[nz][0], n_send, MPI_FLOAT_P, rank + 1, 0, MPI_COMM_WORLD);
        }
        if (has_neighbor_below) {
            MPI_Recv(&array[0][0], n_send, MPI_FLOAT_P, rank - 1, 0, MPI_COMM_WORLD, &status);
        }
    } else {
        if (has_neighbor_below) {
            MPI_Recv(&array[0][0], n_send, MPI_FLOAT_P, rank - 1, 0, MPI_COMM_WORLD, &status);
        }
        if (has_neighbor_above) {
            MPI_Send(&array[nz][0], n_send, MPI_FLOAT_P, rank + 1, 0, MPI_COMM_WORLD);
        }
    }

    // Sending bottom rows to top ghost cells.
    if (rank % 2 == 0) {
        if (has_neighbor_below) {
            MPI_Send(&array[nz_ghost][0], n_send, MPI_FLOAT_P, rank - 1, 0, MPI_COMM_WORLD);
        }
        if (has_neighbor_above) {
            MPI_Recv(&array[nz + nz_ghost][0], n_send, MPI_FLOAT_P, rank + 1, 0, MPI_COMM_WORLD, &status);
        }
    } else {
        if (has_neighbor_above) {
            MPI_Recv(&array[nz + nz_ghost][0], n_send, MPI_FLOAT_P, rank + 1, 0, MPI_COMM_WORLD, &status);
        }
        if (has_neighbor_below) {
            MPI_Send(&array[nz_ghost][0], n_send, MPI_FLOAT_P, rank - 1, 0, MPI_COMM_WORLD);
        }
    }

    // All processes return here if no periodic boundary conditions.
    #if PERIODIC_BOUNDARY_VERTICAL == 0
        return;
    #endif

    if (!mpi_info->has_neighbor_above || !mpi_info->has_neighbor_below)
    {
        communicate_2D_ghost_periodic_boundary(array, mpi_info, nz, nz_ghost, ny);
    }
}