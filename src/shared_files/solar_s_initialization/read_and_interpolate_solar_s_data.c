#include "solar_s_initialization.h"
#include "../array_memory_management/array_memory_management.h"

void read_and_interpolate_solar_s_data(double *r_over_R, double *c_s, double *rho, double *p, double *Gamma_1, double *T, double *H, double *superad_param, double *grad_s, double del_ad, int size)
{
    /*
    Reads the solar_s.h5 file and interpolates the data to the desired grid.

    Parameters
    ----------
    r_over_R : double array
        Pointer to radius array of background state
    c_s : double array
        Pointer to sound speed array of background state
    rho : double array
        Pointer to density array of background state
    p : double array
        Pointer to pressure array of background state
    Gamma_1 : double array
        Pointer to the adiabatic index array of background state
    T : double array
        Pointer to temperature array of background state
    H : double array
        Pointer to pressure scale height array of background state
    superad_param : double array
        Pointer to the superadiabaticity parameter array of background state
    grad_s : double array
        Pointer to the entropy gradient array of background state
    del_ad : double
        Adiabatic temperature gradient
    size : int
        Size of the background state arrays
    */

    int solar_s_size = 2482;
    double *r_over_R_solar_s, *c_s_solar_s, *rho_solar_s, *p_solar_s, *Gamma_1_solar_s, *T_solar_s;
    double *H_solar_s, *superad_param_solar_s, *grad_s_solar_s;

    allocate_1D_array(&r_over_R_solar_s, solar_s_size);
    allocate_1D_array(&c_s_solar_s, solar_s_size);
    allocate_1D_array(&rho_solar_s, solar_s_size);
    allocate_1D_array(&p_solar_s, solar_s_size);
    allocate_1D_array(&Gamma_1_solar_s, solar_s_size);
    allocate_1D_array(&T_solar_s, solar_s_size);
    allocate_1D_array(&H_solar_s, solar_s_size);
    allocate_1D_array(&superad_param_solar_s, solar_s_size);
    allocate_1D_array(&grad_s_solar_s, solar_s_size);

    read_solar_s_data("../additional_files/solar_s.h5", r_over_R_solar_s, c_s_solar_s, rho_solar_s, p_solar_s, Gamma_1_solar_s, T_solar_s, solar_s_size);

    calculate_pressure_scale_height(r_over_R_solar_s, p_solar_s, H_solar_s, solar_s_size);
    calculate_superadiabcicity_parameter(p_solar_s, T_solar_s, superad_param_solar_s, del_ad, solar_s_size);
    calculate_entropy_gradient(p_solar_s, rho_solar_s, T_solar_s, Gamma_1_solar_s, H_solar_s, superad_param_solar_s, grad_s_solar_s, solar_s_size);

    // Interpolate the solar_s data to the grid
    for (int i = 0; i < size; i++)
    {
        interpolate_solar_s(r_over_R_solar_s, c_s_solar_s, rho_solar_s, p_solar_s, Gamma_1_solar_s, T_solar_s, H_solar_s, superad_param_solar_s, grad_s_solar_s, &r_over_R[i], &c_s[i], &rho[i], &p[i], &Gamma_1[i], &T[i], &H[i], &superad_param[i], &grad_s[i], solar_s_size);
    }

    deallocate_1D_array(r_over_R_solar_s);
    deallocate_1D_array(c_s_solar_s);
    deallocate_1D_array(rho_solar_s);
    deallocate_1D_array(p_solar_s);
    deallocate_1D_array(Gamma_1_solar_s);
    deallocate_1D_array(T_solar_s);
    deallocate_1D_array(H_solar_s);
    deallocate_1D_array(superad_param_solar_s);
    deallocate_1D_array(grad_s_solar_s);
}