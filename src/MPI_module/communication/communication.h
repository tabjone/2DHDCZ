#ifndef COMMUNICATION_H__
#define COMMUNICATION_H__

#include "global_float_precision.h"
#include "../mpi_info_struct.h"

void communicate_1D_ghost_above_below(FLOAT_P *array, struct MpiInfo *mpi_info, int nz, int nz_ghost);
void communicate_2D_ghost_above_below(FLOAT_P **array, struct MpiInfo *mpi_info, int nz, int nz_ghost, int ny);
void communicate_3D_ghost_above_below(FLOAT_P ***array, struct MpiInfo *mpi_info, int nz, int nz_ghost, int ny, int nx);

void communicate_1D_ghost_periodic_boundary(FLOAT_P *array, struct MpiInfo *mpi_info, int nz, int nz_ghost);
void communicate_2D_ghost_periodic_boundary(FLOAT_P **array, struct MpiInfo *mpi_info, int nz, int nz_ghost, int ny);
void communicate_3D_ghost_periodic_boundary(FLOAT_P ***array, struct MpiInfo *mpi_info, int nz, int nz_ghost, int ny, int nx);

#endif // COMMUNICATION_H__