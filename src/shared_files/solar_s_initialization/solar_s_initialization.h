#ifndef SOLAR_S_INITIALIZATION_SOLAR_S_INITIALIZATION_H_
#define SOLAR_S_INITIALIZATION_SOLAR_S_INITIALIZATION_H_

#include <hdf5.h>

void read_solar_s_data(const char* filename, double* r_over_R, double* c_s, double* rho, double* p, double* Gamma_1, double* T, hsize_t size);

void interpolate_solar_s(double *r_over_R_solar_s, double *c_s_solar_s, double *rho_solar_s, double *p_solar_s, double *Gamma_1_solar_s, double *T_solar_s, double *r_over_R_i, double *c_s_i, double *rho_i, double *p_i, double *Gamma_1_i, double *T_i, int solar_s_size);

void read_and_interpolate_solar_s_data(double *r_over_R, double *c_s, double *rho, double *p, double *Gamma_1, double *T, int size);

double interpolate_1D(double x0, double x1, double y0, double y1, double x);

#endif // SOLAR_S_INITIALIZATION_SOLAR_S_INITIALIZATION_H_