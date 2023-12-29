#ifndef FIRST_LAW_OF_THERMODYNAMICS_1D_H__
#define FIRST_LAW_OF_THERMODYNAMICS_1D_H__

#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_1D/foreground_variables_struct_1D.h"
#include "data_structures/grid_info/grid_info_1D/grid_info_struct_1D.h"

void first_law_of_thermodynamics_1D(struct ForegroundVariables1D *fg, struct BackgroundVariables *bg, struct GridInfo1D *grid_info);

#endif // FIRST_LAW_OF_THERMODYNAMICS_1D_H__