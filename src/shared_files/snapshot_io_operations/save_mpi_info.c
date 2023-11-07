#include "snapshot_io_operations.h"

void save_mpi_info(struct MpiInfo *mpi_info)
{
    /*
    Saves for now only the total number of processes to a hdf5 file.

    Parameters
    ----------
    mpi_info : struct MpiInfo
        Pointer to the MpiInfo struct.
    */
    
   // Header string
   const char* root_header = "This is the root header";

    // File path
    char file_path[150];
    snprintf(file_path, sizeof(file_path), "%s%s/mpi_info.h5", SAVE_DIR, RUN_NAME);
    #if MPI_ON == 0
    hid_t file;
    herr_t status;

    // Create new default file
    file = H5Fcreate(file_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
        fprintf(stderr, "Failed to create file\n");
    }

    // Adding root header
    add_string_attribute(file, "header", root_header);
    
    // Creating dataspace for scalar
    hid_t dataspace_scalar;

    dataspace_scalar = H5Screate(H5S_SCALAR);
    if (dataspace_scalar < 0) {
        fprintf(stderr, "Failed to create dataspace\n");
    }

    create_write_dataset(file, "total_processes", H5T_NATIVE_INT, dataspace_scalar, &mpi_info->size, "ok");

    // Colose dataspace
    status = H5Sclose(dataspace_scalar);
    if (status < 0) {
        fprintf(stderr, "Failed to close dataspace\n");
    }

    // Closing file
    status = H5Fclose(file);
    if (status < 0) {
        fprintf(stderr, "Failed to close file\n");
    }

    #elif MPI_ON == 1
    if (mpi_info->rank == 0)
    {
    hid_t file;
    herr_t status;

    // Create new default file
    file = H5Fcreate(file_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
        fprintf(stderr, "Failed to create file\n");
    }

    // Adding root header
    add_string_attribute(file, "header", root_header);
    
    // Creating dataspace for scalar
    hid_t dataspace_scalar;

    dataspace_scalar = H5Screate(H5S_SCALAR);
    if (dataspace_scalar < 0) {
        fprintf(stderr, "Failed to create dataspace\n");
    }

    create_write_dataset(file, "total_processes", H5T_NATIVE_INT, dataspace_scalar, &mpi_info->size, "ok");

    // Colose dataspace
    status = H5Sclose(dataspace_scalar);
    if (status < 0) {
        fprintf(stderr, "Failed to close dataspace\n");
    }

    // Closing file
    status = H5Fclose(file);
    if (status < 0) {
        fprintf(stderr, "Failed to close file\n");
    }
    }
    #endif // MPI_ON
}