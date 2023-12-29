#ifndef FIRST_LAW_OF_THERMODYNAMICS_3D_H__
#define FIRST_LAW_OF_THERMODYNAMICS_3D_H__

#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_3D/foreground_variables_struct_3D.h"
#include "data_structures/grid_info/grid_info_3D/grid_info_struct_3D.h"

void first_law_of_thermodynamics_3D(struct ForegroundVariables3D *fg, struct BackgroundVariables *bg, struct GridInfo3D *grid_info);

#endif // FIRST_LAW_OF_THERMODYNAMICS_3D_H__