#include "io_functions.h"
#include "../rhs_functions/rhs_functions.h"

void save_rhs(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables *precalc, int save_nr)
{
    // Create arrays to hold all rhs functions
    FLOAT_P **rhs_vy, **rhs_vz, **rhs_s1;

    // Allocate memory
    allocate_2D_array(&rhs_vy, grid_info->nz, grid_info->ny);
    allocate_2D_array(&rhs_vz, grid_info->nz, grid_info->ny);
    allocate_2D_array(&rhs_s1, grid_info->nz, grid_info->ny);


    // Calculate rhs
    for (int i = 0; i < grid_info->nz; i++)
    {
        for (int j = 0; j < grid_info->ny; j++)
        {
            rhs_vy[i][j] = rhs_dvy_dt_2D(bg, fg, grid_info, precalc, i+grid_info->nz_ghost, j);
            rhs_vz[i][j] = rhs_dvz_dt_2D(bg, fg, grid_info, precalc, i+grid_info->nz_ghost, j);
            rhs_s1[i][j] = rhs_ds1_dt_2D(bg, fg, grid_info, precalc, i+grid_info->nz_ghost, j);
        }
    }

    // File path
    char file_path[150];

    #if MPI_ON == 0
        snprintf(file_path, sizeof(file_path), "%s%s/rhs/rhs%d.h5", SAVE_DIR, RUN_NAME, save_nr);
    #elif MPI_ON == 1
        snprintf(file_path, sizeof(file_path), "%s%s/rhs/rhs%d_%d.h5", SAVE_DIR, RUN_NAME, save_nr, mpi_info->rank);
    #endif // MPI_ON

    // Create file with no groups, but only three datasets 2D
    hsize_t dims[2] = {grid_info->nz, grid_info->ny};
    hid_t file, dataspace_2d;
    herr_t status;

    // Create new default file
    file = H5Fcreate(file_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
        fprintf(stderr, "Failed to create file\n");
    }

    // Create dataspace for 2D arrays
    dataspace_2d = H5Screate_simple(2, dims, NULL);
    if (dataspace_2d < 0) {
        fprintf(stderr, "Failed to create dataspace\n");
    }

    // Create datasets
    create_write_dataset(file, "rhs_vy", H5_FLOAT_P, dataspace_2d, rhs_vy[0], "m/s^2");
    create_write_dataset(file, "rhs_vz", H5_FLOAT_P, dataspace_2d, rhs_vz[0], "m/s^2");
    create_write_dataset(file, "rhs_s1", H5_FLOAT_P, dataspace_2d, rhs_s1[0], "erg/Ks");

    // Close dataspace
    status = H5Sclose(dataspace_2d);
    if (status < 0) {
        fprintf(stderr, "Failed to close dataspace\n");
    }

    // Close file
    status = H5Fclose(file);
    if (status < 0) {
        fprintf(stderr, "Failed to close file\n");
    }

    // Deallocate memory
    deallocate_2D_array(rhs_vy);
    deallocate_2D_array(rhs_vz);
    deallocate_2D_array(rhs_s1);
}