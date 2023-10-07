#ifndef RHS_FUNCTIONS_H__
#define RHS_FUNCTIONS_H__

#include "shared_files.h"
#include "global_parameters.h"


#if DIMENSIONS == 1
// Entropy
FLOAT_P rhs_ds1_dt_1D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i);
// Momentum
FLOAT_P rhs_dvz_dt_1D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i);
// Elliptic equation
FLOAT_P rhs_elliptic_eq_1D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i);
#elif DIMENSIONS == 2
// Entropy
FLOAT_P rhs_ds1_dt_2D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j);
// Momentum
FLOAT_P rhs_dvy_dt_2D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j);
FLOAT_P rhs_dvz_dt_2D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j);
// Elliptic equation
FLOAT_P rhs_elliptic_eq_2D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j);
#elif DIMENSIONS == 3
// Entropy
FLOAT_P rhs_ds1_dt_3D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j, int k);
// Momentum
FLOAT_P rhs_dvx_dt_3D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j, int k);
FLOAT_P rhs_dvy_dt_3D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j, int k);
FLOAT_P rhs_dvz_dt_3D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j, int k);
// Elliptic equation
FLOAT_P rhs_elliptic_eq_3D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j, int k);
#endif // DIMENSIONS

#endif // RHS_FUNCTIONS_H__