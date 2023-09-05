#ifndef SOLAR_S_INITIALIZATION_SOLAR_S_INITIALIZATION_H_
#define SOLAR_S_INITIALIZATION_SOLAR_S_INITIALIZATION_H_

#include <hdf5.h>
#include "../structs/structs.h"

void read_solar_s_data(const char* filename, double* r_over_R, double* rho, double* p, double* T, hsize_t size);

void solar_s_background_initialization(struct BackgroundVariables *bg);

#endif // SOLAR_S_INITIALIZATION_SOLAR_S_INITIALIZATION_H_