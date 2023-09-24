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
    hid_t dataset_p1, dataset_rho1, dataset_T1, dataset_s1, dataset_vx, dataset_vz;
    hid_t dataset_t, dataset_dx, dataset_dz, dataset_nx, dataset_nz, dataset_nz_ghost, dataset_nz_full;
    hid_t dataset_z0, dataset_z1, dataset_x0, dataset_x1;

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

    dataspace_scalar = H5Screate_simple(1, (hsize_t[]){1}, NULL);

    dataset_t = H5Dcreate2(group_grid_data, "t", H5T_NATIVE_DOUBLE, dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_t, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &time);
    dataset_dx = H5Dcreate2(group_grid_data, "dx", H5T_NATIVE_DOUBLE, dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_dx, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dx);
    dataset_dz = H5Dcreate2(group_grid_data, "dz", H5T_NATIVE_DOUBLE, dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_dz, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dz);
    dataset_nx = H5Dcreate2(group_grid_data, "nx", H5T_NATIVE_INT, dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_nx, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nx);
    dataset_nz = H5Dcreate2(group_grid_data, "nz", H5T_NATIVE_INT, dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_nz, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nz);
    dataset_nz_ghost = H5Dcreate2(group_grid_data, "nz_ghost", H5T_NATIVE_INT, dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_nz_ghost, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nz_ghost);
    dataset_nz_full = H5Dcreate2(group_grid_data, "nz_full", H5T_NATIVE_INT, dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_nz_full, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nz_full);
    dataset_z0 = H5Dcreate2(group_grid_data, "z0", H5T_NATIVE_DOUBLE, dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_z0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &z0);
    dataset_z1 = H5Dcreate2(group_grid_data, "z1", H5T_NATIVE_DOUBLE, dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_z1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &z1);
    dataset_x0 = H5Dcreate2(group_grid_data, "x0", H5T_NATIVE_DOUBLE, dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_x0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &x0);
    dataset_x1 = H5Dcreate2(group_grid_data, "x1", H5T_NATIVE_DOUBLE, dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_x1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &x1);

    status = H5Dclose(dataset_t);
    status = H5Dclose(dataset_dx);
    status = H5Dclose(dataset_dz);
    status = H5Dclose(dataset_nx);
    status = H5Dclose(dataset_nz);
    status = H5Dclose(dataset_nz_ghost);
    status = H5Dclose(dataset_nz_full);
    status = H5Dclose(dataset_z0);
    status = H5Dclose(dataset_z1);
    status = H5Dclose(dataset_x0);
    status = H5Dclose(dataset_x1);
    status = H5Sclose(dataspace_scalar);

    status = H5Gclose(group_grid_data);

    if (status < 0) {
        fprintf(stderr, "Failed to close dataspace\n");
    }
    

    // Create group for variables
    group_variables = H5Gcreate2(file, "/variables", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group_variables < 0) {
        fprintf(stderr, "Failed to create variables group\n");
    }

    // print p1 grid
    for (int i = 0; i < nz_full; i++) {
        for (int j = 0; j < nx; j++) {
            printf("%f ", fg->p1[i][j]);
        }
        printf("\n");
    }

    dataspace_2d = H5Screate_simple(2, dims, NULL);
    dataset_T1 = H5Dcreate2(group_variables, "T1", H5T_NATIVE_DOUBLE, dataspace_2d, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_T1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fg->T1);
    dataset_p1 = H5Dcreate2(group_variables, "p1", H5T_NATIVE_DOUBLE, dataspace_2d, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_p1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fg->p1);
    dataset_rho1 = H5Dcreate2(group_variables, "rho1", H5T_NATIVE_DOUBLE, dataspace_2d, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_rho1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fg->rho1);
    dataset_s1 = H5Dcreate2(group_variables, "s1", H5T_NATIVE_DOUBLE, dataspace_2d, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_s1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fg->s1);
    dataset_vx = H5Dcreate2(group_variables, "vx", H5T_NATIVE_DOUBLE, dataspace_2d, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_vx, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fg->vx);
    dataset_vz = H5Dcreate2(group_variables, "vz", H5T_NATIVE_DOUBLE, dataspace_2d, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_vz, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fg->vz);

    status = H5Dclose(dataset_T1);
    status = H5Dclose(dataset_p1);
    status = H5Dclose(dataset_rho1);
    status = H5Dclose(dataset_s1);
    status = H5Dclose(dataset_vx);
    status = H5Dclose(dataset_vz);
    status = H5Gclose(group_variables);

    status = H5Sclose(dataspace_2d);
    if (status < 0) {
        fprintf(stderr, "Failed to close dataspace\n");
    }

    status = H5Fclose(file);
    if (status < 0) {
        fprintf(stderr, "Failed to close file\n");
    }
}