#include "io_functions.h"
#include "../rhs_functions/rhs_functions.h"

void save_rhs(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables *precalc, int save_nr)
{
    /*
    // Create arrays to hold all rhs functions
    FLOAT_P **rhs_vy, **rhs_vz, **rhs_s1, **rhs_p1;
    FLOAT_P **ds1_dy_vy, **ds1_dz_vz;
    FLOAT_P **dvy_dy_vy, **dvy_dz_vz;
    FLOAT_P **dvz_dy_vy, **dvz_dz_vz;

    // Allocate memory
    allocate_2D_array(&rhs_vy, grid_info->nz, grid_info->ny);
    allocate_2D_array(&rhs_vz, grid_info->nz, grid_info->ny);
    allocate_2D_array(&rhs_s1, grid_info->nz, grid_info->ny);
    allocate_2D_array(&rhs_p1, grid_info->nz, grid_info->ny);
    allocate_2D_array(&ds1_dy_vy, grid_info->nz, grid_info->ny);
    allocate_2D_array(&ds1_dz_vz, grid_info->nz, grid_info->ny);
    allocate_2D_array(&dvy_dy_vy, grid_info->nz, grid_info->ny);
    allocate_2D_array(&dvy_dz_vz, grid_info->nz, grid_info->ny);
    allocate_2D_array(&dvz_dy_vy, grid_info->nz, grid_info->ny);
    allocate_2D_array(&dvz_dz_vz, grid_info->nz, grid_info->ny);


    // Calculate rhs
    for (int i = 0; i < grid_info->nz; i++)
    {
        for (int j = 0; j < grid_info->ny; j++)
        {
            rhs_vy[i][j] = rhs_dvy_dt_2D(bg, fg, grid_info, precalc, i+grid_info->nz_ghost, j);
            rhs_vz[i][j] = rhs_dvz_dt_2D(bg, fg, grid_info, precalc, i+grid_info->nz_ghost, j);
            rhs_s1[i][j] = rhs_ds1_dt_2D(bg, fg, grid_info, precalc, i+grid_info->nz_ghost, j);
            rhs_p1[i][j] = rhs_elliptic_eq_2D(bg, fg, grid_info, precalc, i+grid_info->nz_ghost, j);
            ds1_dy_vy[i][j] = upwind_first_derivative_y(fg->s1, fg->vy, precalc, i+grid_info->nz_ghost, j);
            ds1_dz_vz[i][j] = upwind_first_derivative_z(fg->s1, fg->vz, precalc, i+grid_info->nz_ghost, j);
            dvy_dy_vy[i][j] = upwind_first_derivative_y(fg->vy, fg->vy, precalc, i+grid_info->nz_ghost, j);
            dvy_dz_vz[i][j] = upwind_first_derivative_z(fg->vy, fg->vz, precalc, i+grid_info->nz_ghost, j);
            dvz_dy_vy[i][j] = upwind_first_derivative_y(fg->vz, fg->vy, precalc, i+grid_info->nz_ghost, j);
            dvz_dz_vz[i][j] = upwind_first_derivative_z(fg->vz, fg->vz, precalc, i+grid_info->nz_ghost, j);
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
    create_write_dataset(file, "rhs_p1", H5_FLOAT_P, dataspace_2d, rhs_p1[0], "dyn/cm^3");
    create_write_dataset(file, "ds1_dy_vy", H5_FLOAT_P, dataspace_2d, ds1_dy_vy[0], "erg/Kcm");
    create_write_dataset(file, "ds1_dz_vz", H5_FLOAT_P, dataspace_2d, ds1_dz_vz[0], "erg/Kcm");
    create_write_dataset(file, "dvy_dy_vy", H5_FLOAT_P, dataspace_2d, dvy_dy_vy[0], "1/s");
    create_write_dataset(file, "dvy_dz_vz", H5_FLOAT_P, dataspace_2d, dvy_dz_vz[0], "1/s");
    create_write_dataset(file, "dvz_dy_vy", H5_FLOAT_P, dataspace_2d, dvz_dy_vy[0], "1/s");
    create_write_dataset(file, "dvz_dz_vz", H5_FLOAT_P, dataspace_2d, dvz_dz_vz[0], "1/s");

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
    deallocate_2D_array(rhs_p1);
    deallocate_2D_array(ds1_dy_vy);
    deallocate_2D_array(ds1_dz_vz);
    deallocate_2D_array(dvy_dy_vy);
    deallocate_2D_array(dvy_dz_vz);
    deallocate_2D_array(dvz_dy_vy);
    deallocate_2D_array(dvz_dz_vz);
    */
}