#ifndef RHS_FUNCTIONS_3D_H__
#define RHS_FUNCTIONS_3D_H__

#include "global_parameters.h"
#include "shared_files.h"

#include "../spacial_derivatives_3D/spacial_derivatives_3D.h"

FLOAT_P rhs_ds1_dt_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, int i, int j, int k);
FLOAT_P rhs_dvx_dt_3D_vertical_boundary(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, int i, int j, int k);
FLOAT_P rhs_dvx_dt_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, int i, int j, int k);
FLOAT_P rhs_dvy_dt_3D_vertical_boundary(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, int i, int j, int k);
FLOAT_P rhs_dvy_dt_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, int i, int j, int k);
FLOAT_P rhs_dvz_dt_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, int i, int j, int k);
FLOAT_P rhs_elliptic_eq_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, int i, int j, int k);

#endif // RHS_FUNCTIONS_3D_H__