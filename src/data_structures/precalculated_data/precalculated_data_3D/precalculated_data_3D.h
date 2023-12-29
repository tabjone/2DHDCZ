#ifndef PRECALCULATED_DATA_3D_H__
#define PRECALCULATED_DATA_3D_H__

#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/grid_info/grid_info_3D/grid_info_struct_3D.h"
#include "precalculated_data_struct_3D.h"

void allocate_precalculate_data_struct_3D(struct PrecalculatedVariables3D **pv, int nz_full);
void deallocate_precalculated_data_struct_3D(struct PrecalculatedVariables3D **pv);
void initialize_precalculated_data_3D(struct PrecalculatedVariables3D *pv, struct BackgroundVariables *bg, struct GridInfo3D *grid_info);

#endif // PRECALCULATED_DATA_3D_H__