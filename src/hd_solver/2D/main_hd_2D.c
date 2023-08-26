#include "hdf5.h"
#include "../../shared_files/shared_files.h"
#include "./functions.h"
#include "../../shared_files/solar_s_initialization/solar_s_initialization.h"
//#include <math.h>
#include "./structs/structs.h"


// For printing
#include <stdio.h>

int main_hd_2D(int argc, char *argv[])
{
    // Adiabatic temperature gradient, convective instablity when del_ad > del, should be when del_ad > 2/5
    double del_ad = 0.4;

    // Fixed background thermodynamic variables
    //double *r_over_R, *c_s, *p0, *T0, *rho0, *Gamma_1, *H, *superad_param, *grad_s0, *g;

    // Foreground thermodynamic variables
    //double **s1, **rho1, **p1, **T1, **vx, **vz;

    int nz = 200; // Number of grid points in z-direction
    int nx = 100; // Number of grid points in x-direction

    //printf("Hello world!\n");

    struct BackgroundVariables background_variables;
    struct ForegroundVariables foreground_variables;

    solar_s_background_initialization(&background_variables);

    // Saving data to file
    hid_t file_id;
    herr_t status;
    hsize_t dims[1];
    dims[0] = background_variables.nz;
    file_id = H5Fcreate("../data/solar_s_background_3.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    create_and_write_dataset_1D(file_id, "r", dims, background_variables.r);
    create_and_write_dataset_1D(file_id, "T0", dims, background_variables.T0);
    create_and_write_dataset_1D(file_id, "rho0", dims, background_variables.rho0);
    create_and_write_dataset_1D(file_id, "p0", dims, background_variables.p0);
    create_and_write_dataset_1D(file_id, "g", dims, background_variables.g);
    status = H5Fclose(file_id);

    deallocate_background_struct(&background_variables);

    /*
    allocate_background_struct(nz, &background_variables);
    allocate_foreground_struct(nz, nx, &foreground_variables);

    one_time_step(&background_variables, &foreground_variables, nz, nx);

    deallocate_background_struct(&background_variables);
    deallocate_foreground_struct(&foreground_variables);
    */

    /*
    // Allocating memory for 1D arrays
    allocate_1D_array(&r_over_R, nz);
    allocate_1D_array(&T0, nz);
    allocate_1D_array(&rho0, nz);
    allocate_1D_array(&p0, nz);
    allocate_1D_array(&c_s, nz);
    allocate_1D_array(&Gamma_1, nz);
    allocate_1D_array(&H, nz);
    allocate_1D_array(&superad_param, nz);
    allocate_1D_array(&grad_s0, nz);
    allocate_1D_array(&g, nz);

    // Allocating memory for 2D arrays
    allocate_2D_array(&s1, nz, nx);
    allocate_2D_array(&rho1, nz, nx);
    allocate_2D_array(&p1, nz, nx);
    allocate_2D_array(&T1, nz, nx);
    allocate_2D_array(&vx, nz, nx);
    allocate_2D_array(&vz, nz, nx);


    // Creating the radius array and calculating gravity
    double r_over_R_start = 0.10;
    double r_over_R_end = 0.95;

    for (int i=0; i<nz; i++)
    {
        r_over_R[i] = r_over_R_start + 1.0*i/(nz-1) * (r_over_R_end-r_over_R_start);
        g[i] = 1.0/(r_over_R[i]*r_over_R[i]); // temp value
    }

    read_and_interpolate_solar_s_data(r_over_R, c_s, rho0, p0, Gamma_1, T0, H, superad_param, grad_s0, g, del_ad, nz);

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
    create_and_write_dataset_1D(file_id, "g", dims, g);

    status = H5Fclose(file_id);

    double time_to_run = 10.0;
    double dt = 0.1;
    double t = 0.0;
    double dx = 0.1;
    double dz = 0.1;

    double test1, test2, test3;
    int i = 10;
    int j = 10;
    while (t < time_to_run)
    {
        one_time_step_2D_hd(s1, p1, rho1, T1, vx, vz, grad_s0, p0, rho0, T0, g, nz, nx, dx, dz);
        t += dt;
    }

    // Save background thermodynamic variables to file
    // ....

    // Deallocating 1D arrays
    deallocate_1D_array(T0);
    deallocate_1D_array(grad_s0);
    deallocate_1D_array(rho0);
    deallocate_1D_array(p0);
    deallocate_1D_array(c_s);
    deallocate_1D_array(Gamma_1);
    deallocate_1D_array(r_over_R);
    deallocate_1D_array(H);
    deallocate_1D_array(superad_param);
    deallocate_1D_array(g);

    // Deallocating 2D arrays
    deallocate_2D_array(s1);
    deallocate_2D_array(rho1);
    deallocate_2D_array(p1);
    deallocate_2D_array(T1);
    deallocate_2D_array(vx);
    deallocate_2D_array(vz);
    */
    return 0;
}