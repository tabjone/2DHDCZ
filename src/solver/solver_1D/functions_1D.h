#ifndef FUNCTIONS_1D_H__
#define FUNCTIONS_1D_H__

#include "data_structures/data_structures.h"
#include "MPI_module/mpi_info_struct.h"

void initialize_simulation_1D(struct BackgroundVariables **bg, struct ForegroundVariables1D **fg, struct ForegroundVariables1D **fg_previous, struct GridInfo1D **grid_info, struct PrecalculatedVariables1D **precalc, struct MpiInfo *mpi_info, int *save_nr, FLOAT_P *t, FLOAT_P *first_t);
void load_simulation_1D();

#endif // FUNCTIONS_1D_H__