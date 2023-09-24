#include "io_functions.h"

void save_foreground(struct ForegroundVariables2D *fg, struct GridInfo *grid_info, int snap_number, double time)
{
    char file_path[150];

    snprintf(file_path, sizeof(file_path), "data/%s/snap%d.h5", RUN_NAME, snap_number);

    double z0 = grid_info->z0;
    double z1 = grid_info->z1;
    double x0 = grid_info->x0;
    double x1 = grid_info->x1;
    int nz = grid_info->nz;
    int nz_ghost = grid_info->nz_ghost;
    int nz_full = grid_info->nz_full;
    int nx = grid_info->nx;
    double dx = grid_info->dx;
    double dz = grid_info->dz;

    hid_t file, group_grid_data, group_variables;
    hid_t dataspace_scalar, dataspace_2d;

    hsize_t dims[2] = {nz_full, nx};
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

    create_write_dataset(group_grid_data, "t", H5T_NATIVE_DOUBLE, dataspace_scalar, &time);
    create_write_dataset(group_grid_data, "dx", H5T_NATIVE_DOUBLE, dataspace_scalar, &dx);
    create_write_dataset(group_grid_data, "dz", H5T_NATIVE_DOUBLE, dataspace_scalar, &dz);
    create_write_dataset(group_grid_data, "nx", H5T_NATIVE_INT, dataspace_scalar, &nx);
    create_write_dataset(group_grid_data, "nz", H5T_NATIVE_INT, dataspace_scalar, &nz);
    create_write_dataset(group_grid_data, "nz_ghost", H5T_NATIVE_INT, dataspace_scalar, &nz_ghost);
    create_write_dataset(group_grid_data, "nz_full", H5T_NATIVE_INT, dataspace_scalar, &nz_full);
    create_write_dataset(group_grid_data, "z0", H5T_NATIVE_DOUBLE, dataspace_scalar, &z0);
    create_write_dataset(group_grid_data, "z1", H5T_NATIVE_DOUBLE, dataspace_scalar, &z1);
    create_write_dataset(group_grid_data, "x0", H5T_NATIVE_DOUBLE, dataspace_scalar, &x0);
    create_write_dataset(group_grid_data, "x1", H5T_NATIVE_DOUBLE, dataspace_scalar, &x1);

    status = H5Sclose(dataspace_scalar);
    if (status < 0) {
        fprintf(stderr, "Failed to close dataspace\n");
    }

    status = H5Gclose(group_grid_data);
    if (status < 0) {
        fprintf(stderr, "Failed to close dataspace\n");
    }
    
    // Create group for variables
    group_variables = H5Gcreate2(file, "/variables", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group_variables < 0) {
        fprintf(stderr, "Failed to create variables group\n");
    }

    dataspace_2d = H5Screate_simple(2, dims, NULL);
    if (dataspace_2d < 0) {
        fprintf(stderr, "Failed to create dataspace\n");
    }
    
    create_write_dataset(group_variables, "T1", H5T_NATIVE_DOUBLE, dataspace_2d, fg->T1[0]);
    create_write_dataset(group_variables, "p1", H5T_NATIVE_DOUBLE, dataspace_2d, fg->p1[0]);
    create_write_dataset(group_variables, "rho1", H5T_NATIVE_DOUBLE, dataspace_2d, fg->rho1[0]);
    create_write_dataset(group_variables, "s1", H5T_NATIVE_DOUBLE, dataspace_2d, fg->s1[0]);
    create_write_dataset(group_variables, "vx", H5T_NATIVE_DOUBLE, dataspace_2d, fg->vx[0]);
    create_write_dataset(group_variables, "vz", H5T_NATIVE_DOUBLE, dataspace_2d, fg->vz[0]);

    status = H5Gclose(group_variables);
    if (status < 0) {
        fprintf(stderr, "Failed to close group\n");
    }
    status = H5Sclose(dataspace_2d);
    if (status < 0) {
        fprintf(stderr, "Failed to close dataspace\n");
    }

    status = H5Fclose(file);
    if (status < 0) {
        fprintf(stderr, "Failed to close file\n");
    }
}