#include "io_functions.h"
/*
void load_foreground(struct ForegroundVariables2D **fg, const char *file_path) {
    hid_t file = H5Fopen(file_path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        fprintf(stderr, "Failed to open file\n");
        return;
    }

    // Read scalar variables to first allocate the struct
    int nz, nx;
    double dz, dx;
    hid_t dataset;


    
    hid_t file = H5Fopen(file_path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        fprintf(stderr, "Failed to open file\n");
        return;
    }

    // Reading scalar variables
    hid_t dataset = H5Dopen2(file, "/data/nz", H5P_DEFAULT);
    H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(fg->nz));
    H5Dclose(dataset);

    // Repeat above two lines for nz_ghost, nz_full, nx, dx, and dz

    // Reading 2D arrays. You'll need to allocate memory and then read into it.
    hsize_t dims[2];
    hid_t dataspace;

    
    dataset = H5Dopen2(file, "/variables/p1", H5P_DEFAULT);
    dataspace = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    fg->p1 = (double **)malloc(dims[0] * sizeof(double *));
    for (int i = 0; i < dims[0]; i++) {
        fg->p1[i] = (double *)malloc(dims[1] * sizeof(double));
        // You can use a subset reading function to read row by row, or flatten your 2D array and read all at once, then rearrange.
    }
    H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fg->p1[0]);  // Assuming continuous memory, read all at once
    H5Dclose(dataset);

    // Repeat the block above for rho1, T1, s1, vx, and vz

    H5Fclose(file);
}
*/