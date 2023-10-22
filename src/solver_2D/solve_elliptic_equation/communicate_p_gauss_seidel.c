#include "solve_elliptic_equation.h"
#include <mpi.h>

void communicate_p_gauss_seidel(FLOAT_P **array, struct GridInfo *grid_info, struct MpiInfo *mpi_info)
{
    /*
    Communicates the boundary of the below and above processes to our ghost cells

    Parameters
    ----------
    array : FLOAT_P **
        Pressure field
    grid_info : struct GridInfo
        Grid parameters
    mpi_info : struct MpiInfo
        MPI parameters
    */
    
    MPI_Status status;
    
    // Getting grid info
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
            MPI_Send(&array[nz][0], ny, MPI_FLOAT_P, rank + 1, 0, MPI_COMM_WORLD);
        }
        if (has_neighbor_below) {
            MPI_Recv(array[0], ny, MPI_FLOAT_P, rank - 1, 0, MPI_COMM_WORLD, &status);
        }
    } else {
        if (has_neighbor_below) {
            MPI_Recv(array[0], ny, MPI_FLOAT_P, rank - 1, 0, MPI_COMM_WORLD, &status);
        }
        if (has_neighbor_above) {
            MPI_Send(&array[nz][0], ny, MPI_FLOAT_P, rank + 1, 0, MPI_COMM_WORLD);
        }
    }

    // Now sending to below
    if (rank % 2 == 0) {
        if (has_neighbor_below) {
            MPI_Send(&array[1][0], ny, MPI_FLOAT_P, rank - 1, 0, MPI_COMM_WORLD);
        }
        if (has_neighbor_above) {
            MPI_Recv(&array[nz + 1][0], ny, MPI_FLOAT_P, rank + 1, 0, MPI_COMM_WORLD, &status);
        }
    } else {
        if (has_neighbor_above) {
            MPI_Recv(&array[nz + 1][0], ny, MPI_FLOAT_P, rank + 1, 0, MPI_COMM_WORLD, &status);
        }
        if (has_neighbor_below) {
            MPI_Send(&array[1][0], ny, MPI_FLOAT_P, rank - 1, 0, MPI_COMM_WORLD);
        }
    }
}