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
    FLOAT_P *r;
    FLOAT_P *rho0;
    FLOAT_P *p0;
    FLOAT_P *m;
    FLOAT_P *T0;
    FLOAT_P *s0;
    FLOAT_P *grad_s0;
    int N;
    int N_increment;
};

void allocate_integration_variables(struct IntegrationVariables *i_var);
void deallocate_integration_variables(struct IntegrationVariables *i_var);
void set_initial_integration_values(struct IntegrationVariables *i_var, FLOAT_P r_i, FLOAT_P p0_i, FLOAT_P T0_i, FLOAT_P rho0_i, FLOAT_P m_i);

FLOAT_P get_k_value(FLOAT_P r);

void integrate_one_step(struct IntegrationVariables *bg, int i, bool updown);

void read_solar_s_data(const char* filename, FLOAT_P* r_over_R, FLOAT_P* rho0, FLOAT_P* p0, FLOAT_P* T0, hsize_t size);

void solar_s_background_initialization(struct BackgroundVariables *bg, struct MpiInfo *mpi_info, int grid_nz_full, int grid_nz_ghost, FLOAT_P grid_dz, FLOAT_P grid_z0, FLOAT_P grid_z1, FLOAT_P grid_nz);

#endif // SOLAR_S_INITIALIZATION_SOLAR_S_INITIALIZATION_H_