#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <mpi.h>

#include "functions.h"

int main(int argc, char *argv[])
{
    struct MpiInfo *mpi_info;
    mpi_info = malloc(sizeof(struct MpiInfo));

    // Initialize MPI
    #if MPI_ON == 1
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_info->size);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_info->rank);

        // Check if there are neighbors above and below
        if (mpi_info->rank == 0)
        {
            mpi_info->has_neighbor_below = false;
            mpi_info->has_neighbor_above = true;
        }
        else if (mpi_info->rank == mpi_info->size - 1)
        {
            mpi_info->has_neighbor_below = true;
            mpi_info->has_neighbor_above = false;
        }
        else
        {
            mpi_info->has_neighbor_below = true;
            mpi_info->has_neighbor_above = true;
        }
        #else
            mpi_info->has_neighbor_below = false;
            mpi_info->has_neighbor_above = false;
    #endif // MPI_ON

    char dir_name[100];
    struct stat st = {0}; // Initialize the stat structure

    // Use snprintf to create the directory path with the base directory and RUN_NAME
    snprintf(dir_name, sizeof(dir_name), "data/%s", RUN_NAME);

    #if MPI_ON == 1
        if (mpi_info->rank == 0)
        {
            // Check if directory exists
            if (stat(dir_name, &st) == -1) {
                // If directory doesn't exist, create it with permissions set to rwxr-xr-x
                if (mkdir(dir_name, 0755) == -1) {
                    perror("Error creating directory");
                    return 1;
                }
            }
        }
        // Wait for rank 0 to finish creating the directory
        MPI_Barrier(MPI_COMM_WORLD);
    #else
        // Check if directory exists
        if (stat(dir_name, &st) == -1) {
            // If directory doesn't exist, create it with permissions set to rwxr-xr-x
            if (mkdir(dir_name, 0755) == -1) {
                perror("Error creating directory");
                return 1;
            }
        }
    #endif // MPI_ON

    #if DIMENSIONS == 1
        printf("Not implemented yet\n");
    #elif DIMENSIONS == 2
        main_2D(argc, argv, mpi_info);
    #elif DIMENSIONS == 3
        main_3D(argc, argv, mpi_info);
    #endif // DIMENSIONS

    free(mpi_info);

    #if MPI_ON == 1
        MPI_Finalize();
    #endif // MPI_ON
    return 0;
}
