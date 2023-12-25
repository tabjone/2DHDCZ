#include "boundary.h"

void update_vertical_boundary_ghostcells_2D(FLOAT_P **array, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info)
{
    int nz = grid_info->nz;
    int nz_ghost = grid_info->nz_ghost;
    int ny = grid_info->ny;

    

    #if MPI_ON == 1
        // Sending and receiving ghost cells
        communicate_2D_ghost_above_below(array, mpi_info, nz, nz_ghost, ny);
    #endif // MPI_ON
    
    // All non periodic boundary conditions should extrapolate ghost cells
    #if PERIODIC_BOUNDARY_VERTICAL == 0
        int nz_full = grid_info->nz_full;
        // If not periodic boundary extrapolate ghost cells at top and bottom
        if (!mpi_info->has_neighbor_below)
        {
            extrapolate_2D_array_down(array, nz_ghost, ny);
        }
        if (!mpi_info->has_neighbor_above)
        {
            extrapolate_2D_array_up(array, nz_full, nz_ghost, ny);
        }
    #endif // VERTICAL_BOUNDARY_TYPE

    // Handle case of periodic boundary
    #if PERIODIC_BOUNDARY_VERTICAL == 1
        #if MPI_ON == 1
            MPI_Status status;
            int size = mpi_info->size;
            // Send and recieve to/from bottom
            if (!mpi_info->has_neighbor_above)
            {
                MPI_Send(&array[nz][0], nz_ghost*ny, MPI_FLOAT_P, 0, 0, MPI_COMM_WORLD);
                // Recieve from bottom
                MPI_Recv(&array[nz+nz_ghost][0], nz_ghost*ny, MPI_FLOAT_P, 0, 1, MPI_COMM_WORLD, &status);
            }
            // Send and recieve to/from top
            if (!mpi_info->has_neighbor_below)
            {
                // Recieve from top
                MPI_Recv(&array[0][0], nz_ghost*ny, MPI_FLOAT_P, size-1, 0, MPI_COMM_WORLD, &status);
                // Send to top
                MPI_Send(&array[nz_ghost][0], nz_ghost*ny, MPI_FLOAT_P, size-1, 1, MPI_COMM_WORLD);
            }
        #else
            // If not mpi just do periodic boundary function
            periodic_boundary_2D(array, grid_info);
        #endif // MPI_ON
    #endif // VERTICAL_BOUNDARY_TYPE
}