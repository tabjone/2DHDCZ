#include "mpi_functions.h"

void communicate_2D_ghost_above_below(FLOAT_P **array, struct GridInfo *grid_info, struct MpiInfo *mpi_info) 
{
    /*
    Communicates ghost cells above and below the current process for a 2D array.

    Parameters
    ----------
    array : float pointer
        Pointer to the array to be communicated.
    grid_info : struct GridInfo pointer
        Pointer to the struct containing grid info.
    mpi_info : struct MpiInfo pointer
        Pointer to the struct containing mpi info.
    */
    MPI_Status status;
    
    // Getting grid info
    int nz_full = grid_info->nz_full;
    int nz = grid_info->nz;
    int nz_ghost = grid_info->nz_ghost;
    int ny = grid_info->ny;

    // Getting mpi info
    int rank = mpi_info->rank;
    bool has_neighbor_above = mpi_info->has_neighbor_above;
    bool has_neighbor_below = mpi_info->has_neighbor_below;

    // Communicating ghost cells with even-odd method
    if (rank % 2 == 0) {
        if (has_neighbor_above) {
            MPI_Send(&array[nz][0], nz_ghost*ny, MPI_FLOAT_P, rank + 1, 0, MPI_COMM_WORLD);
        }
        if (has_neighbor_below) {
            MPI_Recv(array[0], nz_ghost*ny, MPI_FLOAT_P, rank - 1, 0, MPI_COMM_WORLD, &status);
        }
    } else {
        if (has_neighbor_below) {
            MPI_Recv(array[0], nz_ghost*ny, MPI_FLOAT_P, rank - 1, 0, MPI_COMM_WORLD, &status);
        }
        if (has_neighbor_above) {
            MPI_Send(&array[nz][0], nz_ghost*ny, MPI_FLOAT_P, rank + 1, 0, MPI_COMM_WORLD);
        }
    }

    // Now sending to below
    if (rank % 2 == 0) {
        if (has_neighbor_below) {
            MPI_Send(&array[nz_ghost][0], nz_ghost*ny, MPI_FLOAT_P, rank - 1, 0, MPI_COMM_WORLD);
        }
        if (has_neighbor_above) {
            MPI_Recv(&array[nz + nz_ghost][0], nz_ghost*ny, MPI_FLOAT_P, rank + 1, 0, MPI_COMM_WORLD, &status);
        }
    } else {
        if (has_neighbor_above) {
            MPI_Recv(&array[nz + nz_ghost][0], nz_ghost*ny, MPI_FLOAT_P, rank + 1, 0, MPI_COMM_WORLD, &status);
        }
        if (has_neighbor_below) {
            MPI_Send(&array[nz_ghost][0], nz_ghost*ny, MPI_FLOAT_P, rank - 1, 0, MPI_COMM_WORLD);
        }
    }
}