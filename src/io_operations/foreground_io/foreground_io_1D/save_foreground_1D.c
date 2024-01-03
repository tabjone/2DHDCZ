#include "data_structures/grid_info/grid_info_1D/grid_info_struct_1D.h"
#include "data_structures/foreground_data/foreground_data_1D/foreground_variables_struct_1D.h"
#include "MPI_module/mpi_info_struct.h"
#include "io_operations/general_io/general_io.h"
#include "global_parameters.h"

void save_foreground_1D(struct ForegroundVariables1D *fg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info, int snap_number, FLOAT_P time)
{
    /*
    Saves the foreground variables to a hdf5 file at the current time step.

    Parameters
    ----------
    fg : struct ForegroundVariables1D
        Pointer to the foreground variables struct.
    grid_info : struct GridInfo1D
        Pointer to the grid info struct.
    snap_number : int
        The current snapshot number.
    time : FLOAT_P
        The current time.
    */

    // Header strings
    const char* root_header = "This is the root header";
    const char* grid_data_header = "This is the grid data header";
    const char* variables_header = "This is the variables header";

    // File path
    char file_path[150];

    #if MPI_ON == 0
        snprintf(file_path, sizeof(file_path), "%s%s/snap%d.h5", SAVE_DIR, RUN_NAME, snap_number);
    #elif MPI_ON == 1
        snprintf(file_path, sizeof(file_path), "%s%s/snap%d_%d.h5", SAVE_DIR, RUN_NAME, snap_number, mpi_info->rank);
    #endif // MPI_ON
  
    // Getting grid info
    FLOAT_P z0 = grid_info->z0;
    FLOAT_P z1 = grid_info->z1;
    FLOAT_P y0 = grid_info->y0;
    FLOAT_P y1 = grid_info->y1;
    int nz = grid_info->nz;
    int nz_ghost = grid_info->nz_ghost;
    int nz_full = grid_info->nz_full;
    int ny = grid_info->ny;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P dz = grid_info->dz;


    hid_t file, group_grid_data, group_variables;
    hid_t dataspace_scalar, dataspace_1D;

    hsize_t dims[2] = {nz_full, ny};
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
        create_write_dataset(group_grid_data, "t", H5_FLOAT_P, dataspace_scalar, &time, "s");
        create_write_dataset(group_grid_data, "dz", H5_FLOAT_P, dataspace_scalar, &dz, "cm");
        create_write_dataset(group_grid_data, "nz", H5T_NATIVE_INT, dataspace_scalar, &nz, "Grid points in z");
        create_write_dataset(group_grid_data, "nz_ghost", H5T_NATIVE_INT, dataspace_scalar, &nz_ghost, "Ghost points in z");
        create_write_dataset(group_grid_data, "nz_full", H5T_NATIVE_INT, dataspace_scalar, &nz_full, "Full grid points in z");
        create_write_dataset(group_grid_data, "z0", H5_FLOAT_P, dataspace_scalar, &z0, "cm");
        create_write_dataset(group_grid_data, "z1", H5_FLOAT_P, dataspace_scalar, &z1, "cm");
    #else
        create_write_dataset(group_grid_data, "t", H5_FLOAT_P, dataspace_scalar, &time, "s");
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

    status = H5Gclose(group_grid_data);
    if (status < 0) {
        fprintf(stderr, "Failed to close dataspace\n");
    }
    
    // Create group for variables
    group_variables = H5Gcreate2(file, "/variables", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group_variables < 0) {
        fprintf(stderr, "Failed to create variables group\n");
    }

    // Add variables header
    add_string_attribute(group_variables, "header", variables_header);

    dataspace_1D = H5Screate_simple(2, dims, NULL);
    if (dataspace_1D < 0) {
        fprintf(stderr, "Failed to create dataspace\n");
    }
    
    #if UNITS == 0
        create_write_dataset(group_variables, "T1", H5_FLOAT_P, dataspace_1D, fg->T1[0], "K");
        create_write_dataset(group_variables, "p1", H5_FLOAT_P, dataspace_1D, fg->p1[0], "dyn/cm^2");
        create_write_dataset(group_variables, "rho1", H5_FLOAT_P, dataspace_1D, fg->rho1[0], "g/cm^3");
        create_write_dataset(group_variables, "s1", H5_FLOAT_P, dataspace_1D, fg->s1[0], "erg/K");
        create_write_dataset(group_variables, "vz", H5_FLOAT_P, dataspace_1D, fg->vz[0], "cm/s");
    #else
        create_write_dataset(group_variables, "T1", H5_FLOAT_P, dataspace_1D, fg->T1[0], "K");
        create_write_dataset(group_variables, "p1", H5_FLOAT_P, dataspace_1D, fg->p1[0], "Pa/m^2");
        create_write_dataset(group_variables, "rho1", H5_FLOAT_P, dataspace_1D, fg->rho1[0], "kg/m^3");
        create_write_dataset(group_variables, "s1", H5_FLOAT_P, dataspace_1D, fg->s1[0], "J/K");
        create_write_dataset(group_variables, "vz", H5_FLOAT_P, dataspace_1D, fg->vz[0], "m/s");
    #endif // UNITS == 0

    status = H5Gclose(group_variables);
    if (status < 0) {
        fprintf(stderr, "Failed to close group\n");
    }
    status = H5Sclose(dataspace_1D);
    if (status < 0) {
        fprintf(stderr, "Failed to close dataspace\n");
    }

    status = H5Fclose(file);
    if (status < 0) {
        fprintf(stderr, "Failed to close file\n");
    }
}