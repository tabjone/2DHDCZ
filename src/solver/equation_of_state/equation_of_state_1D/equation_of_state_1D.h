#ifndef EQUATION_OF_STATE_1D_H__
#define EQUATION_OF_STATE_1D_H__

#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_1D/foreground_variables_struct_1D.h"
#include "data_structures/grid_info/grid_info_1D/grid_info_struct_1D.h"

void equation_of_state_1D(struct ForegroundVariables1D *fg, struct BackgroundVariables *bg, struct GridInfo1D *grid_info);

#endif // EQUATION_OF_STATE_1D_H__