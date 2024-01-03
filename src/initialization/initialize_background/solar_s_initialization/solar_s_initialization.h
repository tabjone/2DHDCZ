#ifndef SOLAR_S_INITIALIZATION_H__
#define SOLAR_S_INITIALIZATION_H__

#include "global_float_precision.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "MPI_module/mpi_info_struct.h"

void get_initial_values_from_solar_s(FLOAT_P *p0_initial, FLOAT_P *T0_initial, FLOAT_P *rho0_initial, FLOAT_P *m_initial, FLOAT_P r_integration_start);
FLOAT_P get_k_value(FLOAT_P r);
FLOAT_P get_rho_ideal_gas(FLOAT_P p0, FLOAT_P T0);
void solar_s_initialization(struct BackgroundVariables *bg, struct MpiInfo *mpi_info, int my_nz_full, int my_nz_ghost, FLOAT_P dz, FLOAT_P my_z0);

#endif // SOLAR_S_INITIALIZATION_H__