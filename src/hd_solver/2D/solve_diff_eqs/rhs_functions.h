#ifndef RHS_FUNCTIONS_H__
#define RHS_FUNCTIONS_H__

#include "shared_files.h"
#include "global_parameters.h"

FLOAT_P rhs_ds1_dt_1D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i);
FLOAT_P rhs_ds1_dt_2D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j);
FLOAT_P rhs_ds1_dt_3D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j, int k);

FLOAT_P rhs_dvz_dt_2D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i);

FLOAT_P rhs_dvy_dt_2D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j);
FLOAT_P rhs_dvz_dt_2D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j);

FLOAT_P rhs_dvy_dt_3D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j, int k);
FLOAT_P rhs_dvx_dt_3D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j, int k);
FLOAT_P rhs_dvz_dt_3D(struct BackgroundVariables *bg, struct ForegroundVariables *fg, struct GridInfo *grid_info, int i, int j, int k);

#endif // RHS_FUNCTIONS_H__