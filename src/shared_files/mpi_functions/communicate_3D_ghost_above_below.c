#include "mpi_functions.h"

void communicate_3D_ghost_above_below(float ***array, struct GridInfo *grid_info, struct MpiInfo *mpi_info) 
{
    /*
    Communicates ghost cells above and below the current process for a 3D array.

    Parameters
    ----------
    array : float pointer
        Pointer to the array to be communicated.
    grid_info : struct GridInfo pointer
        Pointer to the struct containing grid info.
    mpi_info : struct MpiInfo pointer
        Pointer to the struct containing mpi info.
    */
   /*
    MPI_Status status;

    // Getting grid info
    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;
    int ny = grid_info->ny;
    int nx = grid_info->nx;

    // Getting mpi info
    int rank = mpi_info->rank;
    bool has_neighbor_above = mpi_info->has_neighbor_above;
    bool has_neighbor_below = mpi_info->has_neighbor_below;

    // Communicating ghost cells with even-odd method
    if (rank % 2 == 0) {
        if (has_neighbor_above) {
            MPI_Send(&array[nz_full - 2*nz_ghost][0][0], nz_ghost*ny*nx, MPI_FLOAT_P, rank + 1, 0, MPI_COMM_WORLD);
        }
        if (has_neighbor_below) {
            MPI_Recv(&array[0][0][0], nz_ghost*ny*nx, MPI_FLOAT_P, rank - 1, 0, MPI_COMM_WORLD, &status);
        }
    } else {
        if (has_neighbor_below) {
            MPI_Recv(&array[0][0][0], nz_ghost*ny*nx, MPI_FLOAT_P, rank - 1, 0, MPI_COMM_WORLD, &status);
        }
        if (has_neighbor_above) {
            MPI_Send(&array[nz_full - 2*nz_ghost][0][0], nz_ghost*ny*nx, MPI_FLOAT_P, rank + 1, 0, MPI_COMM_WORLD);
        }
    }
    */
}