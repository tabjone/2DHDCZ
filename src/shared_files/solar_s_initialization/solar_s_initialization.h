#ifndef SOLAR_S_INITIALIZATION_SOLAR_S_INITIALIZATION_H_
#define SOLAR_S_INITIALIZATION_SOLAR_S_INITIALIZATION_H_

#include <hdf5.h>

void calculate_pressure_scale_height(double *r_over_R, double *p0, double *H, int nz);
void calculate_superadiabcicity_parameter(double *p0, double *T0, double *superadiabacicity_parameter, double del_ad, int nz);
void calculate_entropy_gradient(double *p0, double *rho0, double *T0, double *Gamma_1, double *H, double *superadiabacicity_parameter, double *entropy_gradient, int nz);

void read_solar_s_data(const char* filename, double* r_over_R, double* c_s, double* rho, double* p, double* Gamma_1, double* T, hsize_t size);

void interpolate_solar_s(double *r_over_R_solar_s, double *c_s_solar_s, double *rho_solar_s, double *p_solar_s, double *Gamma_1_solar_s, double *T_solar_s, double *H_solar_s, double *superad_param_solar_s, double *grad_s0_solar_s, double *r_over_R_i, double *c_s_i, double *rho_i, double *p_i, double *Gamma_1_i, double *T_i, double *H_i, double *superad_param_i, double *grad_s0_i, int solar_s_size);

void read_and_interpolate_solar_s_data(double *r_over_R, double *c_s, double *rho, double *p, double *Gamma_1, double *T, double *H, double *superad_param, double *grad_s, double del_ad, int size);

#endif // SOLAR_S_INITIALIZATION_SOLAR_S_INITIALIZATION_H_