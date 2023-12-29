#ifndef EQUATION_OF_STATE_2D_H__
#define EQUATION_OF_STATE_2D_H__

#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"

void equation_of_state_2D(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info);

#endif // EQUATION_OF_STATE_2D_H__