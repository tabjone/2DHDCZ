#include "hdf5.h"
#include <stdlib.h>

void create_and_write_dataset_1D(hid_t file_id, const char* datasetname, hsize_t* dims, double* data) 
{
    /*
    Creates a HDF5 file and writes a 1D dataset to it.

    Parameters
    ----------
    file_id : hid_t
        The file id of the HDF5 file.
    datasetname : const char*
        The name of the saved file.
    dims : hsize_t*
        The dimensions of the dataset.
    data : double*
        The data to be written to the dataset.
    */
    // Create the dataspace for the dataset.
    hid_t dataspace_id, dataset_id;  
    herr_t status;

    // Create the dataset.
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate2(file_id, datasetname, H5T_NATIVE_DOUBLE, dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // Write the dataset.
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    // Close the dataset and the dataspace.
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
}