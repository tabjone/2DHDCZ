#include "io_functions.h"

void save_foreground(struct ForegroundVariables2D *fg, int snap_number, double time)
{
    char file_path[100];
    char info[256]; // To store the formatted info string


    snprintf(file_path, sizeof(file_path), "data/%s/snap%d.h5", RUN_NAME, snap_number);
    snprintf(info, sizeof(info), 
    "t: %.2f\nnx: %d\nnz: %d\nnz_ghost: %d\ndx: %.2f\ndz: %.2f", time, fg->nx, fg->nz, fg->nz_ghost, fg->dx, fg->dz);

    hid_t file_id = 0;
    hsize_t dims_info[1] = {1}; // for the string info
    hid_t dataspace_id, dataset_id, filetype;


    if (file_id < 0) {
        perror("Failed to create HDF5 file");
    }

    herr_t status;
    hsize_t dims[2] = {fg->nz_full, fg->nx};
    file_id = H5Fcreate(file_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Save the info string to the HDF5 file
    dataspace_id = H5Screate_simple(1, dims_info, NULL);
    filetype = H5Tcopy(H5T_C_S1);
    H5Tset_size(filetype, sizeof(info));
    dataset_id = H5Dcreate2(file_id, "info", filetype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, filetype, H5S_ALL, H5S_ALL, H5P_DEFAULT, info);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    status = H5Tclose(filetype);


    create_and_write_dataset_2D(file_id, "p1", dims, fg->p1);
    create_and_write_dataset_2D(file_id, "rho1", dims, fg->rho1);
    create_and_write_dataset_2D(file_id, "T1", dims, fg->T1);
    create_and_write_dataset_2D(file_id, "s1", dims, fg->s1);
    create_and_write_dataset_2D(file_id, "vx", dims, fg->vx);
    create_and_write_dataset_2D(file_id, "vz", dims, fg->vz);

    status = H5Fclose(file_id);
    if (status < 0) fprintf(stderr, "Failed to close file\n");
}