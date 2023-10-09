#ifndef MPI_FUNCTIONS_H__
#define MPI_FUNCTIONS_H__

#include <mpi.h>
#include "shared_files.h"
#include "global_parameters.h"

#if DIMENSIONS == 1
void communicate_1D_ghost_above_below(float *array, struct GridInfo *grid_info, struct MpiInfo *mpi_info);
#elif DIMENSIONS == 2
void communicate_2D_ghost_above_below(float **array, struct GridInfo *grid_info, struct MpiInfo *mpi_info);
#elif DIMENSIONS == 3
void communicate_3D_ghost_above_below(float ***array, struct GridInfo *grid_info, struct MpiInfo *mpi_info);
#endif // DIMENSIONS

#endif // MPI_FUNCTIONS_H__