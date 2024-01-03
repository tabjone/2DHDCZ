#ifndef BACKGROUND_IO_H__
#define BACKGROUND_IO_H__

#include "global_float_precision.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "hdf5.h"
#include "MPI_module/mpi_info_struct.h"

void read_solar_s_data(const char* filename, FLOAT_P* r_over_R, FLOAT_P* rho0, FLOAT_P* p0, FLOAT_P* T0, hsize_t size);

void save_background(struct BackgroundVariables *bg, struct MpiInfo *mpi_info, int nz_full, int nz, int nz_ghost, FLOAT_P dz, FLOAT_P z0, FLOAT_P z1);

#endif // BACKGROUND_IO_H__