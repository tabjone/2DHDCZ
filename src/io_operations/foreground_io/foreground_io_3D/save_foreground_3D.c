#include "data_structures/grid_info/grid_info_3D/grid_info_struct_3D.h"
#include "data_structures/foreground_data/foreground_data_3D/foreground_variables_struct_3D.h"
#include "MPI_module/mpi_info_struct.h"
#include "io_operations/general_io/general_io.h"
#include "global_parameters.h"

void save_foreground_3D(struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info, int snap_number, FLOAT_P time)
{
    /*
    Saves the foreground variables to a hdf5 file at the current time step.

    Parameters
    ----------
    fg : struct ForegroundVariables2D
        Pointer to the foreground variables struct.
    grid_info : struct GridInfo2D
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

    snprintf(file_path, sizeof(file_path), "%s%s/snap%d_%d.h5", SAVE_DIR, RUN_NAME, snap_number, mpi_info->rank);
  
    // Getting grid info
    FLOAT_P z0 = grid_info->z0;
    FLOAT_P z1 = grid_info->z1;
    FLOAT_P y0 = grid_info->y0;
    FLOAT_P y1 = grid_info->y1;
    FLOAT_P x0 = grid_info->x0;
    FLOAT_P x1 = grid_info->x1;
    int nz = grid_info->nz;
    int nz_ghost = grid_info->nz_ghost;
    int nz_full = grid_info->nz_full;
    int ny = grid_info->ny;
    int nx = grid_info->nx;
    FLOAT_P dz = grid_info->dz;
    FLOAT_P dy = grid_info->dy;
    FLOAT_P dx = grid_info->dx;

    hid_t file, group_grid_data, group_variables;
    hid_t dataspace_scalar, dataspace_3d;

    hsize_t dims[3] = {nz_full, ny, nx};
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
        create_write_dataset(group_grid_data, "dx", H5_FLOAT_P, dataspace_scalar, &dx, "cm");
        create_write_dataset(group_grid_data, "dy", H5_FLOAT_P, dataspace_scalar, &dy, "cm");
        create_write_dataset(group_grid_data, "dz", H5_FLOAT_P, dataspace_scalar, &dz, "cm");
        create_write_dataset(group_grid_data, "nx", H5T_NATIVE_INT, dataspace_scalar, &nx, "Grid points in x");
        create_write_dataset(group_grid_data, "ny", H5T_NATIVE_INT, dataspace_scalar, &ny, "Grid points in y");
        create_write_dataset(group_grid_data, "nz", H5T_NATIVE_INT, dataspace_scalar, &nz, "Grid points in z");
        create_write_dataset(group_grid_data, "nz_ghost", H5T_NATIVE_INT, dataspace_scalar, &nz_ghost, "Ghost points in z");
        create_write_dataset(group_grid_data, "nz_full", H5T_NATIVE_INT, dataspace_scalar, &nz_full, "Full grid points in z");
        create_write_dataset(group_grid_data, "z0", H5_FLOAT_P, dataspace_scalar, &z0, "cm");
        create_write_dataset(group_grid_data, "z1", H5_FLOAT_P, dataspace_scalar, &z1, "cm");
        create_write_dataset(group_grid_data, "y0", H5_FLOAT_P, dataspace_scalar, &y0, "cm");
        create_write_dataset(group_grid_data, "y1", H5_FLOAT_P, dataspace_scalar, &y1, "cm");
        create_write_dataset(group_grid_data, "x0", H5_FLOAT_P, dataspace_scalar, &x0, "cm");
        create_write_dataset(group_grid_data, "x1", H5_FLOAT_P, dataspace_scalar, &x1, "cm");
    #else
        create_write_dataset(group_grid_data, "t", H5_FLOAT_P, dataspace_scalar, &time, "s");
        create_write_dataset(group_grid_data, "dx", H5_FLOAT_P, dataspace_scalar, &dx, "m");
        create_write_dataset(group_grid_data, "dy", H5_FLOAT_P, dataspace_scalar, &dy, "m");
        create_write_dataset(group_grid_data, "dz", H5_FLOAT_P, dataspace_scalar, &dz, "m");
        create_write_dataset(group_grid_data, "nx", H5T_NATIVE_INT, dataspace_scalar, &nx, "Grid points in x");
        create_write_dataset(group_grid_data, "ny", H5T_NATIVE_INT, dataspace_scalar, &ny, "Grid points in y");
        create_write_dataset(group_grid_data, "nz", H5T_NATIVE_INT, dataspace_scalar, &nz, "Grid points in z");
        create_write_dataset(group_grid_data, "nz_ghost", H5T_NATIVE_INT, dataspace_scalar, &nz_ghost, "Ghost points in z");
        create_write_dataset(group_grid_data, "nz_full", H5T_NATIVE_INT, dataspace_scalar, &nz_full, "Full grid points in z");
        create_write_dataset(group_grid_data, "z0", H5_FLOAT_P, dataspace_scalar, &z0, "m");
        create_write_dataset(group_grid_data, "z1", H5_FLOAT_P, dataspace_scalar, &z1, "m");
        create_write_dataset(group_grid_data, "y0", H5_FLOAT_P, dataspace_scalar, &y0, "m");
        create_write_dataset(group_grid_data, "y1", H5_FLOAT_P, dataspace_scalar, &y1, "m");
        create_write_dataset(group_grid_data, "x0", H5_FLOAT_P, dataspace_scalar, &x0, "m");
        create_write_dataset(group_grid_data, "x1", H5_FLOAT_P, dataspace_scalar, &x1, "m");
    #endif // UNITS

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

    dataspace_3d = H5Screate_simple(3, dims, NULL);
    if (dataspace_3d < 0) {
        fprintf(stderr, "Failed to create dataspace\n");
    }
    
    #if UNITS == 0
        create_write_dataset(group_variables, "T1", H5_FLOAT_P, dataspace_3d, fg->T1[0][0], "K");
        create_write_dataset(group_variables, "p1", H5_FLOAT_P, dataspace_3d, fg->p1[0][0], "dyn/cm^2");
        create_write_dataset(group_variables, "rho1", H5_FLOAT_P, dataspace_3d, fg->rho1[0][0], "g/cm^3");
        create_write_dataset(group_variables, "s1", H5_FLOAT_P, dataspace_3d, fg->s1[0][0], "erg/K");
        create_write_dataset(group_variables, "vx", H5_FLOAT_P, dataspace_3d, fg->vx[0][0], "cm/s");
        create_write_dataset(group_variables, "vy", H5_FLOAT_P, dataspace_3d, fg->vy[0][0], "cm/s");
        create_write_dataset(group_variables, "vz", H5_FLOAT_P, dataspace_3d, fg->vz[0][0], "cm/s");
    #else
        create_write_dataset(group_variables, "T1", H5_FLOAT_P, dataspace_3d, fg->T1[0][0], "K");
        create_write_dataset(group_variables, "p1", H5_FLOAT_P, dataspace_3d, fg->p1[0][0], "Pa/m^2");
        create_write_dataset(group_variables, "rho1", H5_FLOAT_P, dataspace_3d, fg->rho1[0][0], "kg/m^3");
        create_write_dataset(group_variables, "s1", H5_FLOAT_P, dataspace_3d, fg->s1[0][0], "J/K");
        create_write_dataset(group_variables, "vx", H5_FLOAT_P, dataspace_3d, fg->vx[0][0], "m/s");
        create_write_dataset(group_variables, "vy", H5_FLOAT_P, dataspace_3d, fg->vy[0][0], "m/s");
        create_write_dataset(group_variables, "vz", H5_FLOAT_P, dataspace_3d, fg->vz[0][0], "m/s");
    #endif // UNITS

    status = H5Gclose(group_variables);
    if (status < 0) {
        fprintf(stderr, "Failed to close group\n");
    }
    status = H5Sclose(dataspace_3d);
    if (status < 0) {
        fprintf(stderr, "Failed to close dataspace\n");
    }

    status = H5Fclose(file);
    if (status < 0) {
        fprintf(stderr, "Failed to close file\n");
    }
}