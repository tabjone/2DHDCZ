#include "io_functions.h"

void save_array(FLOAT_P **array, int nz, int ny, const char file_path)
{
    hsize_t dims[2] = {nz, ny};
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
    create_write_dataset(file, "data", H5_FLOAT_P, dataspace_2d, array[0], "");

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
}