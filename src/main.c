#include <stdio.h>
#include <stdlib.h>
#include "hdf5.h"

#include "functions.h"

int main(int argc, char *argv[])
{   
    /*
    2D is in the z-y plane
    1D is in the z direction
    */

    double dx, dy, dz;

    hid_t file_id;
    herr_t status;
    hsize_t dims[3];

    dims[0] = 20;
    dims[1] = 10;
    dims[2] = 15;

    dx = 1.0/(dims[2]-1);
    dy = 1.0/(dims[1]-1);
    dz = 1.0/(dims[0]-1);


    double ***threeD_test;
    allocate_3D_array(&threeD_test, dims[0], dims[1], dims[2]);
    
    double **twoD_test;
    allocate_2D_array(&twoD_test, dims[0], dims[1]);
    
    double *oneD_test;
    allocate_1D_array(&oneD_test, dims[0]);

    for (int i=0; i < dims[0]; i++)
    {
        oneD_test[i] = i;
        for (int j=0; j < dims[1]; j++)
        {
            twoD_test[i][j] = i;
            for (int k=0; k < dims[2]; k++)
            {
                threeD_test[i][j][k] = i;
            }
        }
    }
    file_id = H5Fcreate("../data/snap0.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    create_and_write_dataset_1D(file_id, "1D_data", dims, oneD_test);
    create_and_write_dataset_2D(file_id, "2D_data", dims, twoD_test);
    create_and_write_dataset_3D(file_id, "3D_data", dims, threeD_test);
    
    deallocate_3D_array(threeD_test);
    deallocate_2D_array(twoD_test);
    deallocate_1D_array(oneD_test);
}
