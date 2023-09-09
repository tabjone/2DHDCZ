#include "io_functions.h"

void save_foreground(struct ForegroundVariables2D *fg, int snap_number)
{
    char file_path[100];

    snprintf(file_path, sizeof(file_path), "data/%s/snap%d.h5", RUN_NAME, snap_number);

    hid_t file_id = 0;

    if (file_id < 0) {
        perror("Failed to create HDF5 file");
    }

    herr_t status;
    hsize_t dims[2] = {fg->nz_full, fg->nx};

    file_id = H5Fcreate(file_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    create_and_write_dataset_2D(file_id, "p1", dims, fg->p1);
    create_and_write_dataset_2D(file_id, "rho1", dims, fg->rho1);
    create_and_write_dataset_2D(file_id, "T1", dims, fg->T1);
    create_and_write_dataset_2D(file_id, "s1", dims, fg->s1);
    create_and_write_dataset_2D(file_id, "vx", dims, fg->vx);
    create_and_write_dataset_2D(file_id, "vz", dims, fg->vz);

    status = H5Fclose(file_id);
    if (status < 0) fprintf(stderr, "Failed to close file\n");
}