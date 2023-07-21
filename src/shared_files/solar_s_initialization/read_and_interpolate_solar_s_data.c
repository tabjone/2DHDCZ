#include "solar_s_initialization.h"
#include "../array_memory_management/array_memory_management.h"

void read_and_interpolate_solar_s_data(double *r_over_R, double *c_s, double *rho, double *p, double *Gamma_1, double *T, int size)
{
    int solar_s_size = 2482;
    double *r_over_R_solar_s, *c_s_solar_s, *rho_solar_s, *p_solar_s, *Gamma_1_solar_s, *T_solar_s;

    allocate_1D_array(&r_over_R_solar_s, solar_s_size);
    allocate_1D_array(&c_s_solar_s, solar_s_size);
    allocate_1D_array(&rho_solar_s, solar_s_size);
    allocate_1D_array(&p_solar_s, solar_s_size);
    allocate_1D_array(&Gamma_1_solar_s, solar_s_size);
    allocate_1D_array(&T_solar_s, solar_s_size);

    read_solar_s_data("../additional_files/solar_s.h5", r_over_R_solar_s, c_s_solar_s, rho_solar_s, p_solar_s, Gamma_1_solar_s, T_solar_s, solar_s_size);

    // Interpolate the solar_s data to the grid
    for (int i = 0; i < size; i++)
    {
        interpolate_solar_s(r_over_R_solar_s, c_s_solar_s, rho_solar_s, p_solar_s, Gamma_1_solar_s, T_solar_s, &r_over_R[i], &c_s[i], &rho[i], &p[i], &Gamma_1[i], &T[i], solar_s_size);
    }

    deallocate_1D_array(r_over_R_solar_s);
    deallocate_1D_array(c_s_solar_s);
    deallocate_1D_array(rho_solar_s);
    deallocate_1D_array(p_solar_s);
    deallocate_1D_array(Gamma_1_solar_s);
    deallocate_1D_array(T_solar_s);
}