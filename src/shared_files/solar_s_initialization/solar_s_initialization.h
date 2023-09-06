#ifndef SOLAR_S_INITIALIZATION_SOLAR_S_INITIALIZATION_H_
#define SOLAR_S_INITIALIZATION_SOLAR_S_INITIALIZATION_H_

#include <hdf5.h>
#include "global_parameters.h"
#include "global_constants.h"
#include "shared_files.h"
#include <math.h>
#include <stdlib.h>

void read_solar_s_data(const char* filename, double* r_over_R, double* rho0, double* p0, double* T0, hsize_t size);

void solar_s_background_initialization(struct BackgroundVariables *bg);

#endif // SOLAR_S_INITIALIZATION_SOLAR_S_INITIALIZATION_H_