#ifndef RHS_FUNCTIONS_2D_H__
#define RHS_FUNCTIONS_2D_H__

#include "global_float_precision.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "data_structures/precalculated_data/precalculated_data_2D/precalculated_data_struct_2D.h"

// Entropy
FLOAT_P rhs_ds1_dt_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct PrecalculatedVariables2D *precalc, int i, int j);
// Momentum
FLOAT_P rhs_dvy_dt_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct PrecalculatedVariables2D *precalc, int i, int j);
FLOAT_P rhs_dvz_dt_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct PrecalculatedVariables2D *precalc, int i, int j);
// Elliptic equation
FLOAT_P rhs_elliptic_eq_2D(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct PrecalculatedVariables2D *precalc, int i, int j);

void save_derivatives_2D(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct PrecalculatedVariables2D *precalc, struct MpiInfo *mpi_info, int snap_number);

#endif // RHS_FUNCTIONS_2D_H__