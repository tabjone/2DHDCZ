#ifndef FIRST_LAW_OF_THERMODYNAMICS_2D_H_
#define FIRST_LAW_OF_THERMODYNAMICS_2D_H_

#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"

void first_law_thermodynamics_2D(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info);

#endif // FIRST_LAW_OF_THERMODYNAMICS_2D_H_