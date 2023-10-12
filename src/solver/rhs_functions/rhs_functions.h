#ifndef RHS_FUNCTIONS_H__
#define RHS_FUNCTIONS_H__

#include "shared_files.h"
#include "global_parameters.h"
#include "math.h"
#include "global_constants.h"

// Entropy
FLOAT_P rhs_ds1_dt_2D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j);
// Momentum
FLOAT_P rhs_dvy_dt_2D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j);
FLOAT_P rhs_dvy_dt_2D_vertical_boundary(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j);
FLOAT_P rhs_dvz_dt_2D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j);
// Elliptic equation
FLOAT_P rhs_elliptic_eq_2D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j);

#endif // RHS_FUNCTIONS_H__