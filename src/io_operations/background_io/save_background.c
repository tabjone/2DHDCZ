#include "global_float_precision.h"
#include <stdio.h>
#include "data_structures/background_data/background_variables_struct.h"
#include "MPI_module/mpi_info_struct.h"
#include "global_parameters.h"
#include "../general_io/general_io.h"

void save_background(struct BackgroundVariables *bg, struct MpiInfo *mpi_info, int nz_full, int nz, int nz_ghost, FLOAT_P dz, FLOAT_P z0, FLOAT_P z1)
{
    const char* root_header = "This is the root header";
    const char* grid_data_header = "This is the grid data header";
    const char* variables_header = "This is the variables header";

    // This should be re-written as foreground saving
    char file_path[150];
    snprintf(file_path, sizeof(file_path), "%s%s/background_%d.h5", SAVE_DIR, RUN_NAME, mpi_info->rank);

    hid_t file, group_grid_data, group_variables;
    hid_t dataspace_scalar, dataspace_1d;

    hsize_t dims[1] = {nz_full};
    herr_t status;

    // Create new default file
    file = H5Fcreate(file_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
        fprintf(stderr, "Failed to create file\n");
    }

    // Adding root header
    add_string_attribute(file, "header", root_header);

    // Create group for grid_data
    group_grid_data = H5Gcreate2(file, "/grid_info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group_grid_data < 0) {
        fprintf(stderr, "Failed to create grid_data group\n");
    }

    // Add grid_data header
    add_string_attribute(group_grid_data, "header", grid_data_header);

    dataspace_scalar = H5Screate(H5S_SCALAR);
    if (dataspace_scalar < 0) {
        fprintf(stderr, "Failed to create dataspace\n");
    }

    #if UNITS == 0
        create_write_dataset(group_grid_data, "dz", H5_FLOAT_P, dataspace_scalar, &dz, "cm");
        create_write_dataset(group_grid_data, "nz", H5T_NATIVE_INT, dataspace_scalar, &nz, "Grid points in z");
        create_write_dataset(group_grid_data, "nz_ghost", H5T_NATIVE_INT, dataspace_scalar, &nz_ghost, "Ghost points in z");
        create_write_dataset(group_grid_data, "nz_full", H5T_NATIVE_INT, dataspace_scalar, &nz_full, "Full grid points in z");
        create_write_dataset(group_grid_data, "z0", H5_FLOAT_P, dataspace_scalar, &z0, "cm");
        create_write_dataset(group_grid_data, "z1", H5_FLOAT_P, dataspace_scalar, &z1, "cm");
    #else
        create_write_dataset(group_grid_data, "dz", H5_FLOAT_P, dataspace_scalar, &dz, "m");
        create_write_dataset(group_grid_data, "nz", H5T_NATIVE_INT, dataspace_scalar, &nz, "Grid points in z");
        create_write_dataset(group_grid_data, "nz_ghost", H5T_NATIVE_INT, dataspace_scalar, &nz_ghost, "Ghost points in z");
        create_write_dataset(group_grid_data, "nz_full", H5T_NATIVE_INT, dataspace_scalar, &nz_full, "Full grid points in z");
        create_write_dataset(group_grid_data, "z0", H5_FLOAT_P, dataspace_scalar, &z0, "m");
        create_write_dataset(group_grid_data, "z1", H5_FLOAT_P, dataspace_scalar, &z1, "m");
    #endif // UNITS == 0


    status = H5Sclose(dataspace_scalar);
    if (status < 0) {
        fprintf(stderr, "Failed to close dataspace\n");
    }

    // Close group
    status = H5Gclose(group_grid_data);
    if (status < 0) {
        fprintf(stderr, "Failed to close group\n");
    }

    // Create group for variables
    group_variables = H5Gcreate2(file, "/variables", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Add variables header
    add_string_attribute(group_variables, "header", variables_header);

    // Create 1D dataspace
    dataspace_1d = H5Screate_simple(1, dims, NULL);

    #if UNITS == 0
        create_write_dataset(group_variables, "r", H5_FLOAT_P, dataspace_1d, bg->r, "cm");
        create_write_dataset(group_variables, "T0", H5_FLOAT_P, dataspace_1d, bg->T0, "K");
        create_write_dataset(group_variables, "rho0", H5_FLOAT_P, dataspace_1d, bg->rho0, "g/cm^3");
        create_write_dataset(group_variables, "p0", H5_FLOAT_P, dataspace_1d, bg->p0, "dyn/cm^2");
        create_write_dataset(group_variables, "g", H5_FLOAT_P, dataspace_1d, bg->g, "cm/s^2");
        create_write_dataset(group_variables, "grad_s0", H5_FLOAT_P, dataspace_1d, bg->grad_s0, "erg/(cm K)");
    #else
        create_write_dataset(group_variables, "r", H5_FLOAT_P, dataspace_1d, bg->r, "m");
        create_write_dataset(group_variables, "T0", H5_FLOAT_P, dataspace_1d, bg->T0, "K");
        create_write_dataset(group_variables, "rho0", H5_FLOAT_P, dataspace_1d, bg->rho0, "kg/m^3");
        create_write_dataset(group_variables, "p0", H5_FLOAT_P, dataspace_1d, bg->p0, "Pa/m^2");
        create_write_dataset(group_variables, "g", H5_FLOAT_P, dataspace_1d, bg->g, "m/s^2");
        create_write_dataset(group_variables, "grad_s0", H5_FLOAT_P, dataspace_1d, bg->grad_s0, "J/K");
    #endif // UNITS == 0


    status = H5Sclose(dataspace_1d);
    status = H5Gclose(group_variables);

    // Close file
    status = H5Fclose(file);
    if (status < 0) {
        fprintf(stderr, "Failed to close file\n");
    }
}