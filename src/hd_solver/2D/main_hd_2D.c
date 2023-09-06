#include "hdf5.h"
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>

#include "functions.h"


int main_hd_2D(int argc, char *argv[])
{
    int nx = x_size/dx - 1; // Number of grid points in x-direction
    int nz = (R_END - R_START)/dz - 1; // Number of grid points in z-direction 

    struct BackgroundVariables background_variables;

    struct ForegroundVariables2D *foreground_variables, *foreground_variables_previous, *tmp_ptr;
    foreground_variables = (struct ForegroundVariables2D *)malloc(sizeof(struct ForegroundVariables2D));
    foreground_variables_previous = (struct ForegroundVariables2D *)malloc(sizeof(struct ForegroundVariables2D));

    allocate_background_struct(nz, &background_variables);
    allocate_foreground_struct_2D(nz, nx, foreground_variables);
    allocate_foreground_struct_2D(nz, nx, foreground_variables_previous);

    solar_s_background_initialization(&background_variables);

    char file_path[150];
    int k = 0;

    // Construct the full path for the snapshot file inside the new directory
    snprintf(file_path, sizeof(file_path), "data/%s/snap%d.h5", RUN_NAME, k);

    // Saving data to file
    hid_t file_id = 0;

    if (file_id < 0) {
        perror("Failed to create HDF5 file");
    }
    
    herr_t status;
    hsize_t dims[1];
    dims[0] = background_variables.nz_full;
    file_id = H5Fcreate(file_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    create_and_write_dataset_1D(file_id, "r", dims, background_variables.r);
    create_and_write_dataset_1D(file_id, "T0", dims, background_variables.T0);
    create_and_write_dataset_1D(file_id, "rho0", dims, background_variables.rho0);
    create_and_write_dataset_1D(file_id, "p0", dims, background_variables.p0);
    create_and_write_dataset_1D(file_id, "g", dims, background_variables.g);
    create_and_write_dataset_1D(file_id, "grad_s0", dims, background_variables.grad_s0);
    status = H5Fclose(file_id);
    if (status < 0) fprintf(stderr, "Failed to close file\n");
    
    initialize_foreground_struct_zeros(foreground_variables_previous);
    // Copy foreground_variables to foreground_temp
    deep_copy_foreground_2D(foreground_variables, foreground_variables_previous);


    //compute foreground_variables using foreground_variables_previous
    one_time_step(&background_variables, foreground_variables_previous, foreground_variables);
    // pointer swap
    tmp_ptr = foreground_variables_previous;
    foreground_variables_previous = foreground_variables;
    foreground_variables = tmp_ptr;

    deallocate_background_struct(&background_variables);
    deallocate_foreground_struct_2D(foreground_variables);
    deallocate_foreground_struct_2D(foreground_variables_previous);
    free(foreground_variables);
    free(foreground_variables_previous);
    return 0;
}