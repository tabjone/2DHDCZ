#ifndef PRECALCULATED_DATA_1D_H__
#define PRECALCULATED_DATA_1D_H__

#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/grid_info/grid_info_1D/grid_info_struct_1D.h"
#include "precalculated_data_struct_1D.h"

void allocate_precalculate_data_struct_1D(struct PrecalculatedVariables1D **pv, int nz_full);
void deallocate_precalculated_data_struct_1D(struct PrecalculatedVariables1D **pv);
void initialize_precalculated_data_1D(struct PrecalculatedVariables1D *pv, struct BackgroundVariables *bg, struct GridInfo1D *grid_info);

#endif // PRECALCULATED_DATA_1D_H__