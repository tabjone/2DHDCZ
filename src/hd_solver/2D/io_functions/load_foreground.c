#include "io_functions.h"

double load_foreground(struct ForegroundVariables2D *fg, struct GridInfo *grid_info, const char *file_path) {
    double time;  
    herr_t status;

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

    hid_t dataset_rho1 = H5Dopen(fg_group, "rho1", H5P_DEFAULT);
    // Read the datasets into the structure
    for (int i = 0; i < grid_info->nz_full; i++)
    {
        if (H5Dread(dataset_rho1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fg->rho1[i]) < 0)
            fprintf(stderr, "Failed to read row %d of rho1 dataset\n", i);
    }
    H5Dclose(dataset_rho1);

    hid_t dataset_T1 = H5Dopen(fg_group, "T1", H5P_DEFAULT);
    for (int i = 0; i < grid_info->nz_full; i++)
    {
        if (H5Dread(dataset_T1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fg->T1[i]) < 0)
            fprintf(stderr, "Failed to read row %d of T1 dataset\n", i);
    }
    H5Dclose(dataset_T1);

    hid_t dataset_p1 = H5Dopen(fg_group, "p1", H5P_DEFAULT);
    for (int i = 0; i < grid_info->nz_full; i++)
    {
        if (H5Dread(dataset_p1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fg->p1[i]) < 0)
            fprintf(stderr, "Failed to read row %d of p1 dataset\n", i);
    }
    H5Dclose(dataset_p1);

    hid_t dataset_s1 = H5Dopen(fg_group, "s1", H5P_DEFAULT);
    for (int i = 0; i < grid_info->nz_full; i++)
    {
        if (H5Dread(dataset_s1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fg->s1[i]) < 0)
            fprintf(stderr, "Failed to read row %d of s1 dataset\n", i);
    }
    H5Dclose(dataset_s1);

    hid_t dataset_vx = H5Dopen(fg_group, "vx", H5P_DEFAULT);
    for (int i = 0; i < grid_info->nz_full; i++)
    {
        if (H5Dread(dataset_vx, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fg->vx[i]) < 0)
            fprintf(stderr, "Failed to read row %d of vx dataset\n", i);
    }
    H5Dclose(dataset_vx);

    hid_t dataset_vz = H5Dopen(fg_group, "vz", H5P_DEFAULT);
    for (int i = 0; i < grid_info->nz_full; i++)
    {
        if (H5Dread(dataset_vz, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fg->vz[i]) < 0)
            fprintf(stderr, "Failed to read row %d of vz dataset\n", i);
    }
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
    status = H5Dread(H5Dopen(grid_info_group, "t", H5P_DEFAULT), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &time);
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
