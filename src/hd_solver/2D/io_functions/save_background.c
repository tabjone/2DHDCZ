#include "io_functions.h"

void save_background(struct BackgroundVariables *bg)
{
    char file_path[150];

    // Construct the full path for the snapshot file inside the new directory
    snprintf(file_path, sizeof(file_path), "data/%s/background.h5", RUN_NAME);

    // Saving data to file
    hid_t file_id = 0;

    if (file_id < 0) {
        perror("Failed to create HDF5 file");
    }
    
    herr_t status;
    hsize_t dims[1] = {bg->nz_full};
    
    file_id = H5Fcreate(file_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    create_and_write_dataset_1D(file_id, "r", dims, bg->r);
    create_and_write_dataset_1D(file_id, "T0", dims, bg->T0);
    create_and_write_dataset_1D(file_id, "rho0", dims, bg->rho0);
    create_and_write_dataset_1D(file_id, "p0", dims, bg->p0);
    create_and_write_dataset_1D(file_id, "g", dims, bg->g);
    create_and_write_dataset_1D(file_id, "grad_s0", dims, bg->grad_s0);
    status = H5Fclose(file_id);
    if (status < 0) fprintf(stderr, "Failed to close file\n");
}