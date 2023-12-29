#ifndef RHS_FUNCTIONS_3D_H__
#define RHS_FUNCTIONS_3D_H__

#include "global_float_precision.h"
#include "spacial_derivatives_module/derivatives_3D/derivatives_3D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_3D/foreground_variables_struct_3D.h"
#include "data_structures/grid_info/grid_info_3D/grid_info_struct_3D.h"
#include "data_structures/precalculated_data/precalculated_data_3D/precalculated_data_3D.h"

// Entropy
FLOAT_P rhs_ds1_dt_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct PrecalculatedVariables3D *precalc, int i, int j, int k);
// Momentum
FLOAT_P rhs_dvx_dt_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct PrecalculatedVariables3D *precalc, int i, int j, int k);
FLOAT_P rhs_dvy_dt_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct PrecalculatedVariables3D *precalc, int i, int j, int k);
FLOAT_P rhs_dvz_dt_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct PrecalculatedVariables3D *precalc, int i, int j, int k);
// Elliptic equation
FLOAT_P rhs_elliptic_eq_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct PrecalculatedVariables3D *precalc, int i, int j, int k);

#endif // RHS_FUNCTIONS_3D_H__