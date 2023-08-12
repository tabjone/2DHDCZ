#include "hdf5.h"
#include "../../shared_files/shared_files.h"
#include "./hd_2D_functions.h"
#include <math.h>

int main_hd_2D(int argc, char *argv[])
{
    // Testing time-step
    one_time_step_2D_hd(10, 10);


    // Adiabatic temperature gradient, convective instablity when del_ad > del, should be when del_ad > 2/5
    double del_ad = 0.4;

    // Fixed background thermodynamic variables
    double *r_over_R, *c_s, *p0, *T0, *rho0, *Gamma_1, *H, *superad_param, *grad_s0;

    int nz = 50; // Number of grid points in z-direction

    // Allocating memory
    allocate_1D_array(&r_over_R, nz);
    allocate_1D_array(&T0, nz);
    allocate_1D_array(&rho0, nz);
    allocate_1D_array(&p0, nz);
    allocate_1D_array(&c_s, nz);
    allocate_1D_array(&Gamma_1, nz);
    allocate_1D_array(&H, nz);
    allocate_1D_array(&superad_param, nz);
    allocate_1D_array(&grad_s0, nz);

    double r_over_R_start = 0.10;
    double r_over_R_end = 0.95;

    for (int i=0; i<nz; i++)
    {
        r_over_R[i] = r_over_R_start + 1.0*i/(nz-1) * (r_over_R_end-r_over_R_start);
    }

    read_and_interpolate_solar_s_data(r_over_R, c_s, rho0, p0, Gamma_1, T0, H, superad_param, grad_s0, del_ad, nz);

    /*
    printf("Read and interpolated solar_s data\n");
    for (int i = 0; i < nz; i++)
    {
        printf("r_over_R=%.2f, c_s=%.2e, rho=%.2e, p=%.2e, Gamma_1=%.2f, T=%.2e\n", r_over_R[i], c_s[i], rho0[i], p0[i], Gamma_1[i], T0[i]);
    }
    */
   
    // Saving data to file
    hid_t file_id;
    herr_t status;
    hsize_t dims[1];
    dims[0] = nz;
    file_id = H5Fcreate("../data/solar_s_background_2.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


    create_and_write_dataset_1D(file_id, "H", dims, H);
    create_and_write_dataset_1D(file_id, "superadiabacicity_parameter", dims, superad_param);
    create_and_write_dataset_1D(file_id, "grad_s0", dims, grad_s0);
    create_and_write_dataset_1D(file_id, "r_over_R", dims, r_over_R);
    create_and_write_dataset_1D(file_id, "c_s", dims, c_s);
    create_and_write_dataset_1D(file_id, "Gamma_1", dims, Gamma_1);
    create_and_write_dataset_1D(file_id, "T0", dims, T0);
    create_and_write_dataset_1D(file_id, "rho0", dims, rho0);
    create_and_write_dataset_1D(file_id, "p0", dims, p0);

    status = H5Fclose(file_id);

    deallocate_1D_array(H);
    deallocate_1D_array(superad_param);



    double time_to_run = 10.0;
    double dt = 0.1;
    double t = 0.0;

    while (t < time_to_run)
    {
        

        // Fixed dt for now
        t += dt;
    }

    // Save background thermodynamic variables to file
    // ....

    deallocate_1D_array(T0);
    deallocate_1D_array(grad_s0);
    deallocate_1D_array(rho0);
    deallocate_1D_array(p0);
    deallocate_1D_array(c_s);
    deallocate_1D_array(Gamma_1);
    deallocate_1D_array(r_over_R);

    return 0;
}