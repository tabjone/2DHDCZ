#ifndef FUNCTIONS_3D_H__
#define FUNCTIONS_3D_H__

#include "data_structures/data_structures.h"
#include "MPI_module/mpi_info_struct.h"

void initialize_simulation_3D(struct BackgroundVariables **bg, struct ForegroundVariables3D **fg, struct ForegroundVariables3D **fg_previous, struct GridInfo3D **grid_info, struct PrecalculatedVariables3D **precalc, struct MpiInfo *mpi_info, int *save_nr, FLOAT_P *t, FLOAT_P *first_t);
void load_simulation_3D();

#endif // FUNCTIONS_3D_H__