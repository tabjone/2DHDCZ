#ifndef FUNCTIONS_2D_H__
#define FUNCTIONS_2D_H__

#include "data_structures/data_structures.h"
#include "MPI_module/mpi_info_struct.h"

void initialize_simulation_2D(struct BackgroundVariables **bg, struct ForegroundVariables2D **fg, struct ForegroundVariables2D **fg_previous, struct GridInfo2D **grid_info, struct PrecalculatedVariables2D **precalc, struct MpiInfo *mpi_info, int *save_nr, FLOAT_P *t, FLOAT_P *first_t);
void load_simulation_2D();

#endif // FUNCTIONS_2D_H__