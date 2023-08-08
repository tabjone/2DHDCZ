#include "hdf5.h"
#include "../../shared_files/shared_files.h"
#include <math.h>

int main_hd_2D(int argc, char *argv[])
{
    // Adiabatic temperature gradient, convective instablity when del_ad > del, should be when del_ad > 2/5
    double del_ad = 0.4;

    // Fixed background thermodynamic variables
    double *r_over_R, *s0, *c_s, *p0, *T0, *rho0, *Gamma_1;

    int nz = 100; // Number of grid points in z-direction

    // Allocating memory
    allocate_1D_array(&r_over_R, nz);
    allocate_1D_array(&s0, nz);
    allocate_1D_array(&T0, nz);
    allocate_1D_array(&rho0, nz);
    allocate_1D_array(&p0, nz);
    allocate_1D_array(&c_s, nz);
    allocate_1D_array(&Gamma_1, nz);

    for (int i=0; i<nz; i++)
    {
        r_over_R[i] = i * 1.0f/(nz-1);
    }

    read_and_interpolate_solar_s_data(r_over_R, c_s, rho0, p0, Gamma_1, T0, nz);
    /*
    printf("Read and interpolated solar_s data\n");
    for (int i = 0; i < nz; i++)
    {
        printf("r_over_R=%.2f, c_s=%.2e, rho=%.2e, p=%.2e, Gamma_1=%.2f, T=%.2e\n", r_over_R[i], c_s[i], rho0[i], p0[i], Gamma_1[i], T0[i]);
    }
    */

    // RHS of eq. (72) in Lantz & Fan (1999)
    double *RHS_background_entropy;
    // Spesific heat with constant pressure
    double c_p;
    // Pressure scale height
    double *H;
    // Superadiabacicity parameter
    double *superadiabacicity_parameter;
    // Ideal gas constant normalized by the average mass per particle
    double r_star;
    // Radius of the sun in cgs units
    double R_sun = 6.957e10;

    // First we calculate the pressure scale height, by definition right after eq. (7) in Lantz & Fan (1999)
    allocate_1D_array(&H, nz);

    // Handle the end points using forward and backward difference
    H[0] = - R_sun * (r_over_R[1] - r_over_R[0]) * p0[0]/ (p0[1] - p0[0]);
    H[nz-1] = - R_sun * (r_over_R[nz-1] - r_over_R[nz-2]) * p0[nz-1]/ (p0[nz-1] - p0[nz-2]);

    // Loop over the rest of the points using central difference
    for (int i = 1; i < nz-1; ++i) {
        H[i] = - R_sun * (r_over_R[i+1] - r_over_R[i-1]) * p0[i]/ (p0[i+1] - p0[i-1]);
    }

    // Then we calculate the superadiabacicity parameter
    allocate_1D_array(&superadiabacicity_parameter, nz);

    // dlnP = dP/P, dlnT = dT/T
    // Handle the end points using forward and backward difference
    superadiabacicity_parameter[0] = - del_ad + p0[0]/T0[0] * (T0[1] - T0[0]) / (p0[1] - p0[0]);
    superadiabacicity_parameter[nz-1] = - del_ad + p0[nz-1]/T0[nz-1] * (T0[nz-1] - T0[nz-2]) / (p0[nz-1] - p0[nz-2]);
    
    // Loop over the rest of the points using central difference
    for (int i = 1; i < nz-1; ++i) {
        superadiabacicity_parameter[i] = - del_ad + p0[i]/T0[i] * (T0[i+1] - T0[i-1]) / (p0[i+1] - p0[i-1]);
    }

    // Then we calculate the RHS of eq. (72) in Lantz & Fan (1999)
    allocate_1D_array(&RHS_background_entropy, nz);

    for (int i = 0; i < nz; i++)
    {
        r_star = p0[i] / rho0[i] / T0[i];
        c_p = r_star / (1 - 1/Gamma_1[i]);

        RHS_background_entropy[i] = - (c_p / H[i]) * superadiabacicity_parameter[i];
    }

    // Integrating eq. (72) in Lantz & Fan (1999) w.r.t z to get the background entropy

    // End points
    r_star = p0[0] / rho0[0] / T0[0];
    s0[0] = r_star / (1 - 1/Gamma_1[0]);
    r_star = p0[nz-1] / rho0[nz-1] / T0[nz-1];
    s0[nz-1] = r_star / (1 - 1/Gamma_1[nz-1]);

    // Loop over the rest of the points using Euler-Forward
    double dz;
    for (int i = 1; i < nz-1; ++i) {
        dz = R_sun * (r_over_R[i] - r_over_R[i-1]);
        s0[i] = s0[i-1] + RHS_background_entropy[i] * dz;
    }


    // Saving data to file
    hid_t file_id;
    herr_t status;
    hsize_t dims[1];
    dims[0] = nz;
    file_id = H5Fcreate("../data/solar_s_background.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    create_and_write_dataset_1D(file_id, "H", dims, H);
    create_and_write_dataset_1D(file_id, "superadiabacicity_parameter", dims, superadiabacicity_parameter);
    create_and_write_dataset_1D(file_id, "s0", dims, s0);
    create_and_write_dataset_1D(file_id, "r_over_R", dims, r_over_R);
    create_and_write_dataset_1D(file_id, "c_s", dims, c_s);
    create_and_write_dataset_1D(file_id, "Gamma_1", dims, Gamma_1);
    create_and_write_dataset_1D(file_id, "T0", dims, T0);
    create_and_write_dataset_1D(file_id, "rho0", dims, rho0);
    create_and_write_dataset_1D(file_id, "p0", dims, p0);


    deallocate_1D_array(H);
    deallocate_1D_array(superadiabacicity_parameter);
    deallocate_1D_array(RHS_background_entropy);


    // Save background thermodynamic variables to file
    // ....

    deallocate_1D_array(T0);
    deallocate_1D_array(s0);
    deallocate_1D_array(rho0);
    deallocate_1D_array(p0);
    deallocate_1D_array(c_s);
    deallocate_1D_array(Gamma_1);
    deallocate_1D_array(r_over_R);

    return 0;
}