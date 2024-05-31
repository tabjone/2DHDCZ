#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <mpi.h>
#include <hdf5.h>
#include "global_parameters.h"
#include <stdlib.h>
#include "functions.h"

#include "MPI_module/MPI_module.h"

int main_run(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    //H5open();
    //hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    //H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    
    struct MpiInfo *mpi_info;
    initialize_mpi_info_struct(&mpi_info);

    char dir_name[100];
    struct stat st = {0}; // Initialize the stat structure

    // Use snprintf to create the directory path with the base directory and RUN_NAME
    snprintf(dir_name, sizeof(dir_name), "%s%s", SAVE_DIR, RUN_NAME);
    
    if (mpi_info->rank == 0)
    {
        // Check if directory exists
        if (stat(dir_name, &st) == -1) {
            // If directory doesn't exist, create it with permissions set to rwxr-xr-x
            if (mkdir(dir_name, 0755) == -1) {
                printf("Error creating directory\n");
                return 1;
            }
        }
    }

    #if SAVE_RHS == 1
        char dir_name_rhs[100];
        char dir_name_elliptic[100];

        // Use snprintf to create the directory path with the base directory and RUN_NAME
        snprintf(dir_name_rhs, sizeof(dir_name_rhs), "%s%s/rhs", SAVE_DIR, RUN_NAME);
        snprintf(dir_name_elliptic, sizeof(dir_name_elliptic), "%s%s/elliptic_vars", SAVE_DIR, RUN_NAME);
        
        // Create directory
        if (mpi_info->rank == 0)
        {
            // Check if directory exists
            if (stat(dir_name_rhs, &st) == -1) {
                // If directory doesn't exist, create it with permissions set to rwxr-xr-x
                if (mkdir(dir_name_rhs, 0755) == -1) {
                    perror("Error creating directory");
                    return 1;
                }
            }

            // Check if directory exists
            if (stat(dir_name_elliptic, &st) == -1) {
                // If directory doesn't exist, create it with permissions set to rwxr-xr-x
                if (mkdir(dir_name_elliptic, 0755) == -1) {
                    perror("Error creating directory");
                    return 1;
                }
            }
        }
    #endif // SAVE_RHS

    // Wait for rank 0 to finish creating the directory
    MPI_Barrier(MPI_COMM_WORLD);
    
    #if DIMENSIONS == 1
        main_1D(argc, argv, mpi_info);
    #elif DIMENSIONS == 2
        main_2D(argc, argv, mpi_info);
    #elif DIMENSIONS == 3
        main_3D(argc, argv, mpi_info);
    #endif // DIMENSIONS
	
	MPI_Barrier(MPI_COMM_WORLD);
    
    free(mpi_info);

    //H5Pclose(plist_id);
    //H5close();
    MPI_Finalize();
    
    return 0;
}
