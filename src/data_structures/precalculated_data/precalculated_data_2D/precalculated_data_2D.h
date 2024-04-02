#ifndef PRECALCULATED_DATA_2D_H__
#define PRECALCULATED_DATA_2D_H__

#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "precalculated_data_struct_2D.h"
#include "MPI_module/mpi_info_struct.h"

void allocate_precalculate_data_struct_2D(struct PrecalculatedVariables2D **pv, int nz_full);
void deallocate_precalculated_data_struct_2D(struct PrecalculatedVariables2D **pv);
void initialize_precalculated_data_2D(struct PrecalculatedVariables2D *pv, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info);

#endif // PRECALCULATED_DATA_2D_H__