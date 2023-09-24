#include "io_functions.h"

void save_background(struct BackgroundVariables *bg, struct GridInfo *grid_info)
{
    // This should be re-written as foreground saving
    char file_path[150];

    // Construct the full path for the snapshot file inside the new directory
    snprintf(file_path, sizeof(file_path), "data/%s/background.h5", RUN_NAME);

    int nz_full = grid_info->nz_full;
    int nz_ghost = grid_info->nz_ghost;
    int nz = grid_info->nz;
    double dz = grid_info->dz;

    hid_t file, group_grid_data, group_variables;
    hid_t dataspace_scalar, dataspace_1d;

    hsize_t dims[1] = {grid_info->nz_full};
    herr_t status;

    // Create new default file
    file = H5Fcreate(file_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
        fprintf(stderr, "Failed to create file\n");
    }

    // Create group for grid_data
    group_grid_data = H5Gcreate2(file, "/grid_info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group_grid_data < 0) {
        fprintf(stderr, "Failed to create grid_data group\n");
    }

    dataspace_scalar = H5Screate(H5S_SCALAR);
    if (dataspace_scalar < 0) {
        fprintf(stderr, "Failed to create dataspace\n");
    }

    create_write_dataset(group_grid_data, "dz", H5T_NATIVE_DOUBLE, dataspace_scalar, &dz);
    create_write_dataset(group_grid_data, "nz", H5T_NATIVE_INT, dataspace_scalar, &nz);
    create_write_dataset(group_grid_data, "nz_ghost", H5T_NATIVE_INT, dataspace_scalar, &nz_ghost);
    create_write_dataset(group_grid_data, "nz_full", H5T_NATIVE_INT, dataspace_scalar, &nz_full);
    create_write_dataset(group_grid_data, "z0", H5T_NATIVE_DOUBLE, dataspace_scalar, &(grid_info->z0));
    create_write_dataset(group_grid_data, "z1", H5T_NATIVE_DOUBLE, dataspace_scalar, &(grid_info->z1));

    status = H5Sclose(dataspace_scalar);
    if (status < 0) {
        fprintf(stderr, "Failed to close dataspace\n");
    }

    // Close group
    status = H5Gclose(group_grid_data);
    if (status < 0) {
        fprintf(stderr, "Failed to close group\n");
    }

    // Create group for variables
    group_variables = H5Gcreate2(file, "/variables", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Create 1D dataspace
    dataspace_1d = H5Screate_simple(1, dims, NULL);

    create_write_dataset(group_variables, "r", H5T_NATIVE_DOUBLE, dataspace_1d, bg->r[0]);
    create_write_dataset(group_variables, "T0", H5T_NATIVE_DOUBLE, dataspace_1d, bg->T0[0]);
    create_write_dataset(group_variables, "rho0", H5T_NATIVE_DOUBLE, dataspace_1d, bg->rho0[0]);
    create_write_dataset(group_variables, "p0", H5T_NATIVE_DOUBLE, dataspace_1d, bg->p0[0]);
    create_write_dataset(group_variables, "g", H5T_NATIVE_DOUBLE, dataspace_1d, bg->g[0]);
    create_write_dataset(group_variables, "grad_s0", H5T_NATIVE_DOUBLE, dataspace_1d, bg->grad_s0[0]);

    status = H5Sclose(dataspace_1d);
    status = H5Gclose(group_variables);

    // Close file
    status = H5Fclose(file);
    if (status < 0) {
        fprintf(stderr, "Failed to close file\n");
    }
}