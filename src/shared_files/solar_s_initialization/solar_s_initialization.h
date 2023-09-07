#ifndef SOLAR_S_INITIALIZATION_SOLAR_S_INITIALIZATION_H_
#define SOLAR_S_INITIALIZATION_SOLAR_S_INITIALIZATION_H_

#include <hdf5.h>
#include "global_parameters.h"
#include "global_constants.h"
#include "shared_files.h"
#include <math.h>
#include <stdlib.h>

struct IntegrationVariables
{
    double *r;
    double *rho0;
    double *p0;
    double *m;
    double *T0;
    double *s0;
    double *grad_s0;
    int N;
    int N_increment;
};

void allocate_integration_variables(struct IntegrationVariables *i_var);
void deallocate_integration_variables(struct IntegrationVariables *i_var);
void set_initial_integration_values(struct IntegrationVariables *i_var, double r_i, double p0_i, double T0_i, double rho0_i, double m_i);

double get_k_value(double r);

void integrate_one_step(struct IntegrationVariables *bg, int i, bool updown);

void read_solar_s_data(const char* filename, double* r_over_R, double* rho0, double* p0, double* T0, hsize_t size);

void solar_s_background_initialization(struct BackgroundVariables *bg);

#endif // SOLAR_S_INITIALIZATION_SOLAR_S_INITIALIZATION_H_