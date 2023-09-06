#include "hdf5.h"
#include <stdlib.h>

void read_dataset_2D(hid_t file_id, const char* datasetname, hsize_t* dims, double** data) 
{
    /*
    Reads and writes a 2D HDF5 dataset into a 2D array.

    Parameters
    ----------
    file_id : hid_t
        File identifier.
    datasetname : const char*
        Name of the dataset to read.
    dims : hsize_t*
        Array of dimensions of the dataset.
    data : double**
        Pointer to the 2D array to write the dataset into.
    */

    // Open the dataset
    hid_t dataset_id = H5Dopen2(file_id, datasetname, H5P_DEFAULT);
    // Read the dataset into the array
    herr_t status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(data)[0][0]);
    // Close the dataset
    status = H5Dclose(dataset_id);
    if (status < 0) fprintf(stderr, "Failed to close file\n");
}