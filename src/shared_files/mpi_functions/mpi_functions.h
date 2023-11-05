#ifndef MPI_FUNCTIONS_H__
#define MPI_FUNCTIONS_H__

#include <mpi.h>
#include "shared_files.h"
#include "global_parameters.h"

void communicate_1D_ghost_above_below(FLOAT_P *array, struct MpiInfo *mpi_info, int nz_full, int nz, int nz_ghost); void communicate_2D_ghost_above_below(FLOAT_P **array, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info);
void communicate_3D_ghost_above_below(FLOAT_P ***array, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info);

void communicate_background_ghost_above_below(struct BackgroundVariables *bg, struct MpiInfo *mpi_info, int nz_full, int nz, int nz_ghost);
#endif // MPI_FUNCTIONS_H__