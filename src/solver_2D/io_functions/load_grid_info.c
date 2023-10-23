#include "io_functions.h"

void load_grid_info(struct GridInfo2D **grid_info, const char *file_path)
{
    /*
    HELLO FUTURE ME, I HAVE PLACED A HUGE BUG HERE BY SETTING THE Z_OFFSET TO 0.0. THIS SHOULD BE SAVED AND LOADED FROM THE FILE.
    */

    // Numbers to load into the grid_info struct
    int nz, nz_ghost, nz_full, ny;
    FLOAT_P dz, dy, z0, z1, y0, y1;

    FLOAT_P z_offset = 0.0;

    herr_t status;

    // Open the file
    hid_t file = H5Fopen(file_path, H5F_ACC_RDONLY, H5P_DEFAULT);

    
    // Read datasets from the grid_info group
    hid_t grid_group = H5Gopen(file, "grid_info", H5P_DEFAULT);
    if (grid_group < 0) {
        fprintf(stderr, "Failed to open grid_info group\n");
        H5Fclose(file);
    }

    H5Dread(H5Dopen(grid_group, "ny", H5P_DEFAULT), H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ny);
    H5Dread(H5Dopen(grid_group, "nz", H5P_DEFAULT), H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nz);
    H5Dread(H5Dopen(grid_group, "nz_ghost", H5P_DEFAULT), H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nz_ghost);
    H5Dread(H5Dopen(grid_group, "nz_full", H5P_DEFAULT), H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nz_full);
    H5Dread(H5Dopen(grid_group, "dz", H5P_DEFAULT), H5_FLOAT_P, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dz);
    H5Dread(H5Dopen(grid_group, "dy", H5P_DEFAULT), H5_FLOAT_P, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dy);
    H5Dread(H5Dopen(grid_group, "z0", H5P_DEFAULT), H5_FLOAT_P, H5S_ALL, H5S_ALL, H5P_DEFAULT, &z0);
    H5Dread(H5Dopen(grid_group, "z1", H5P_DEFAULT), H5_FLOAT_P, H5S_ALL, H5S_ALL, H5P_DEFAULT, &z1);
    H5Dread(H5Dopen(grid_group, "y0", H5P_DEFAULT), H5_FLOAT_P, H5S_ALL, H5S_ALL, H5P_DEFAULT, &y0);
    H5Dread(H5Dopen(grid_group, "y1", H5P_DEFAULT), H5_FLOAT_P, H5S_ALL, H5S_ALL, H5P_DEFAULT, &y1);

    // Close the grid_info group
    status = H5Gclose(grid_group);
    if (status < 0) {
        fprintf(stderr, "Failed to close grid_info group\n");
    }

    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wincompatible-pointer-types"
    // Load the numbers into the grid_info struct
    allocate_grid_info_struct_2D(grid_info, nz, nz_ghost, nz_full, ny, dz, dy, z0, z1, z_offset, y0, y1);
    #pragma GCC diagnostic pop
    // Close the file
    status = H5Fclose(file);
    if (status < 0) {
        fprintf(stderr, "Failed to close file\n");
    }
}