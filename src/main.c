#include <stdio.h>
#include <stdlib.h>
#include "hdf5.h"

#include "functions.h"

int main(int argc, char *argv[])
{
    int size = 10;
    double *r_over_R, *c_s, *rho, *p, *Gamma_1, *T;
    allocate_1D_array(&r_over_R, size);
    allocate_1D_array(&c_s, size);
    allocate_1D_array(&rho, size);
    allocate_1D_array(&p, size);
    allocate_1D_array(&Gamma_1, size);
    allocate_1D_array(&T, size);

    for (int i=0; i<size; i++)
    {
        r_over_R[i] = 0.1*i/(size-1);
    }

    read_and_interpolate_solar_s_data(r_over_R, c_s, rho, p, Gamma_1, T, size);
    printf("Read and interpolated solar_s data\n");
    for (int i = 0; i < size; i++)
    {
        printf("r_over_R=%f, c_s=%f, rho=%f, p=%f, Gamma_1=%f, T=%f\n", r_over_R[i], c_s[i], rho[i], p[i], Gamma_1[i], T[i]);
    }

    deallocate_1D_array(r_over_R);
    deallocate_1D_array(c_s);
    deallocate_1D_array(rho);
    deallocate_1D_array(p);
    deallocate_1D_array(Gamma_1);
    deallocate_1D_array(T);


    /*
    int solar_s_data_size = 2482;
    double *T, *rho, *p, *c_s, *r_over_R, *Gamma_1;
    allocate_1D_array(&T, solar_s_data_size);
    allocate_1D_array(&rho, solar_s_data_size);
    allocate_1D_array(&p, solar_s_data_size);
    allocate_1D_array(&c_s, solar_s_data_size);
    allocate_1D_array(&r_over_R, solar_s_data_size);
    allocate_1D_array(&Gamma_1, solar_s_data_size);

    //read_solar_s_data("../additional_files/solar_s.h5", r_over_R, c_s, rho, p, Gamma_1, T);
    
    const char *filename = "../additional_files/solar_s.h5";
    read_solar_s_data(filename, r_over_R, c_s, rho, p, Gamma_1, T, solar_s_data_size);

    for (int i=0; i < 10; i++)
    {
        printf("%f %f %.10f %f %f %f\n", r_over_R[i], c_s[i], rho[i], p[i], Gamma_1[i], T[i]);
    }

    deallocate_1D_array(T);
    deallocate_1D_array(rho);
    deallocate_1D_array(p);
    deallocate_1D_array(c_s);
    deallocate_1D_array(r_over_R);
    deallocate_1D_array(Gamma_1);

    */
    /*
    2D is in the z-y plane
    1D is in the z direction
    */
   /*
    printf("%f\n", central_first_derivative_second_order(1,2,1));

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
    deallocate_1D_array(oneD_test);*/
}
