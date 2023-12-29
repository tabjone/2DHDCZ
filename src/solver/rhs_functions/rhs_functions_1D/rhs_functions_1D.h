#ifndef RHS_FUNCTIONS_1D_H__
#define RHS_FUNCTIONS_1D_H__

#include "global_float_precision.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_1D/foreground_variables_struct_1D.h"
#include "data_structures/grid_info/grid_info_1D/grid_info_struct_1D.h"
#include "data_structures/precalculated_data/precalculated_data_1D/precalculated_data_struct_1D.h"

FLOAT_P rhs_ds1_dt_1D(struct BackgroundVariables *bg, struct ForegroundVariables1D *fg, struct GridInfo1D *grid_info, struct PrecalculatedVariables1D *precalc, int i);
FLOAT_P rhs_dvz_dt_1D(struct BackgroundVariables *bg, struct ForegroundVariables1D *fg, struct GridInfo1D *grid_info, struct PrecalculatedVariables1D *precalc, int i);
FLOAT_P rhs_elliptic_eq_1D(struct BackgroundVariables *bg, struct ForegroundVariables1D *fg, struct GridInfo1D *grid_info, struct PrecalculatedVariables1D *precalc, int i);

#endif // RHS_FUNCTIONS_1D_H__