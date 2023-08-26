#include "hdf5.h"

void read_solar_s_data(const char* filename, double* r_over_R, double* c_s, double* rho, double* p, double* Gamma_1, double* T, hsize_t size)
{
    /*
    Read the solar s data from a hdf5 file and store it in the arrays passed as arguments.

    Parameters
    ----------
    filename : const char*
        The name of the file to read from.
    r_over_R : double*
        The array to store the r_over_R data in (Normalized solar radius). In units of solar radii.
    c_s : double*
        The array to store the c_s data in (sound speed). In units of cm/s.
    rho : double*
        The array to store the rho data in (density). In units of g/cm^3.
    p : double*
        The array to store the p data in (pressure). In units of dyn/cm^2.
    Gamma_1 : double*
        The array to store the Gamma_1 data in.
    T : double*
        The array to store the T data in (temperature). In units of K.
    */

    hid_t file_id;
    hid_t dataset_id;
    herr_t status;

    /* Open an existing file. */
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "Failed to open file: %s\n", filename);
        return;
    }

    /* Read datasets. */
    dataset_id = H5Dopen2(file_id, "/r_over_R", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, r_over_R);
    H5Dclose(dataset_id);

    dataset_id = H5Dopen2(file_id, "/c_s", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, c_s);
    H5Dclose(dataset_id);

    dataset_id = H5Dopen2(file_id, "/rho", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rho);
    H5Dclose(dataset_id);

    dataset_id = H5Dopen2(file_id, "/p", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, p);
    H5Dclose(dataset_id);

    dataset_id = H5Dopen2(file_id, "/Gamma_1", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Gamma_1);
    H5Dclose(dataset_id);

    dataset_id = H5Dopen2(file_id, "/T", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, T);
    H5Dclose(dataset_id);

    /* Close the file. */
    status = H5Fclose(file_id);
    if (status < 0) fprintf(stderr, "Failed to close file\n");
}