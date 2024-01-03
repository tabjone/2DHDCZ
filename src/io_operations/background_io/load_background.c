#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "hdf5.h"

void load_background(struct BackgroundVariables *bg, struct GridInfo2D *grid_info, const char *file_path) 
{
    // Open the file
    hid_t file = H5Fopen(file_path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        fprintf(stderr, "Failed to open file\n");
        return;
    }

    // Read datasets from the background group
    hid_t bg_group = H5Gopen(file, "variables", H5P_DEFAULT);
    if (bg_group < 0) {
        fprintf(stderr, "Failed to open background group\n");
        H5Fclose(file);
    }

    hid_t dataset_rho0 = H5Dopen(bg_group, "rho0", H5P_DEFAULT);
    if (H5Dread(dataset_rho0, H5_FLOAT_P, H5S_ALL, H5S_ALL, H5P_DEFAULT, bg->rho0) < 0)
        fprintf(stderr, "Failed to read rho0 dataset\n");
    H5Dclose(dataset_rho0);

    hid_t dataset_T0 = H5Dopen(bg_group, "T0", H5P_DEFAULT);
    if (H5Dread(dataset_T0, H5_FLOAT_P, H5S_ALL, H5S_ALL, H5P_DEFAULT, bg->T0) < 0)
        fprintf(stderr, "Failed to read T0 dataset\n");
    H5Dclose(dataset_T0);

    hid_t dataset_p0 = H5Dopen(bg_group, "p0", H5P_DEFAULT);
    if (H5Dread(dataset_p0, H5_FLOAT_P, H5S_ALL, H5S_ALL, H5P_DEFAULT, bg->p0) < 0)
        fprintf(stderr, "Failed to read p0 dataset\n");
    H5Dclose(dataset_p0);

    hid_t dataset_r = H5Dopen(bg_group, "r", H5P_DEFAULT);
    if (H5Dread(dataset_r, H5_FLOAT_P, H5S_ALL, H5S_ALL, H5P_DEFAULT, bg->r) < 0)
        fprintf(stderr, "Failed to read r dataset\n");
    H5Dclose(dataset_r);

    hid_t dataset_grad_s0 = H5Dopen(bg_group, "grad_s0", H5P_DEFAULT);
    if (H5Dread(dataset_grad_s0, H5_FLOAT_P, H5S_ALL, H5S_ALL, H5P_DEFAULT, bg->grad_s0) < 0)
        fprintf(stderr, "Failed to read grad_s0 dataset\n");
    H5Dclose(dataset_grad_s0);

    hid_t dataset_g = H5Dopen(bg_group, "g", H5P_DEFAULT);
    if (H5Dread(dataset_g, H5_FLOAT_P, H5S_ALL, H5S_ALL, H5P_DEFAULT, bg->g) < 0)
        fprintf(stderr, "Failed to read g dataset\n");
    H5Dclose(dataset_g);

    // Close the background group
    H5Gclose(bg_group);
    if (bg_group < 0) {
        fprintf(stderr, "Failed to close background group\n");
    }

    // Close the file
    H5Fclose(file);
    if (file < 0) {
        fprintf(stderr, "Failed to close file\n");
        return;
    }
}