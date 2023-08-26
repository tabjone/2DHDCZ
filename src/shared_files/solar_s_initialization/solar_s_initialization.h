#ifndef SOLAR_S_INITIALIZATION_SOLAR_S_INITIALIZATION_H_
#define SOLAR_S_INITIALIZATION_SOLAR_S_INITIALIZATION_H_

#include <hdf5.h>
#include "../../hd_solver/2D/structs/structs.h"

void calculate_pressure_scale_height(double *r_over_R, double *p0, double *H, int nz);
void calculate_superadiabcicity_parameter(double *p0, double *T0, double *superadiabacicity_parameter, double del_ad, int nz);
void calculate_entropy_gradient(double *p0, double *rho0, double *T0, double *Gamma_1, double *H, double *superadiabacicity_parameter, double *entropy_gradient, int nz);
void calculate_gravitational_acceleration(double *r_over_R, double *rho, double *g, int nz);

void read_solar_s_data(const char* filename, double* r_over_R, double* c_s, double* rho, double* p, double* Gamma_1, double* T, hsize_t size);

void solar_s_background_initialization(struct BackgroundVariables *bg);

#endif // SOLAR_S_INITIALIZATION_SOLAR_S_INITIALIZATION_H_