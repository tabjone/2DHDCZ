#ifndef SOLVE_DIFF_EQS_H__
#define SOLVE_DIFF_EQS_H__

#include "shared_files.h"
#include "global_parameters.h"

double rhs_dvx_dt(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, int i, int j);

double rhs_dvz_dt(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, int i, int j);

double rhs_ds1_dt(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, int i, int j);

//double rhs_ds1_dt_vertical_boundary(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, int i, int j);

double rhs_dvx_dt_horizontal_boundary(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, int i, int j);

#endif // SOLVE_DIFF_EQS_H__