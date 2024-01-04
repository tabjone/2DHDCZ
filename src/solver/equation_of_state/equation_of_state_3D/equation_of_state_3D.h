#ifndef EQUATION_OF_STATE_3D_H__
#define EQUATION_OF_STATE_3D_H__

#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_3D/foreground_variables_struct_3D.h"
#include "data_structures/grid_info/grid_info_3D/grid_info_struct_3D.h"

void equation_of_state_3D(struct ForegroundVariables3D *fg, struct BackgroundVariables *bg, struct GridInfo3D *grid_info);

#endif // EQUATION_OF_STATE_3D_H__