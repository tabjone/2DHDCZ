#include "io_functions.h"
#include "global_parameters.h"

FLOAT_P load_foreground(struct ForegroundVariables2D *fg, struct GridInfo *grid_info, const char *file_path) 
{    
    // There needs to be a function to read a dataset into a 2D array for a given address

    FLOAT_P time;  
    herr_t status;
    hid_t dataset_rho1, dataset_p1, dataset_s1, dataset_vx, dataset_vz;

    // Open the file
    hid_t file = H5Fopen(file_path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        fprintf(stderr, "Failed to open file\n");
        exit(1);
    }

    // Open the variables group
    hid_t fg_group = H5Gopen(file, "variables", H5P_DEFAULT);
    if (fg_group < 0) {
        fprintf(stderr, "Failed to open variables group\n");
        H5Fclose(file);
        exit(1);
    }

    hid_t dataset_T1 = H5Dopen(fg_group, "T1", H5P_DEFAULT);
    if (H5Dread(dataset_T1, H5_FLOAT_P, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fg->T1[0][0]) < 0)
    fprintf(stderr, "Failed to read the T1 dataset\n");

    H5Dclose(dataset_T1);
    
    dataset_rho1 = H5Dopen(fg_group, "rho1", H5P_DEFAULT);
    if (H5Dread(dataset_rho1, H5_FLOAT_P, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fg->rho1[0][0]) < 0)
    fprintf(stderr, "Failed to read the rho1 dataset\n");

    H5Dclose(dataset_rho1);

    dataset_p1 = H5Dopen(fg_group, "p1", H5P_DEFAULT);
    if (H5Dread(dataset_p1, H5_FLOAT_P, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fg->p1[0][0]) < 0)
    fprintf(stderr, "Failed to read the p1 dataset\n");

    H5Dclose(dataset_p1);

    dataset_s1 = H5Dopen(fg_group, "s1", H5P_DEFAULT);
    if (H5Dread(dataset_s1, H5_FLOAT_P, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fg->s1[0][0]) < 0)
    fprintf(stderr, "Failed to read the s1 dataset\n");

    H5Dclose(dataset_s1);

    dataset_vx = H5Dopen(fg_group, "vx", H5P_DEFAULT);
    if (H5Dread(dataset_vx, H5_FLOAT_P, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fg->vx[0][0]) < 0)
    fprintf(stderr, "Failed to read the vx dataset\n");

    H5Dclose(dataset_vx);

    dataset_vz = H5Dopen(fg_group, "vz", H5P_DEFAULT);
    if (H5Dread(dataset_vz, H5_FLOAT_P, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fg->vz[0][0]) < 0)
    fprintf(stderr, "Failed to read the vz dataset\n");

    H5Dclose(dataset_vz);
 
    // Close the variables group
    H5Gclose(fg_group);

    // Open the grid_info group
    hid_t grid_info_group = H5Gopen(file, "grid_info", H5P_DEFAULT);
    if (grid_info_group < 0) {
        fprintf(stderr, "Failed to open grid_info group\n");
        H5Fclose(file);
        exit(1);
    }

    // Read the time
    status = H5Dread(H5Dopen(grid_info_group, "t", H5P_DEFAULT), H5_FLOAT_P, H5S_ALL, H5S_ALL, H5P_DEFAULT, &time);
    if (status < 0) {
        fprintf(stderr, "Failed to read time from grid_info group\n");
        H5Gclose(grid_info_group);
        H5Fclose(file);
        exit(1);
    }

    // Close the grid_info group
    H5Gclose(grid_info_group);

    // Close the file
    H5Fclose(file);

    return time;
}
