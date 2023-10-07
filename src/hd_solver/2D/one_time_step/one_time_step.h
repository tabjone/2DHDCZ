#ifndef ONE_TIME_STEP_H__
#define ONE_TIME_STEP_H__

#include "../functions.h"
#include "shared_files.h"

FLOAT_P rk1(struct BackgroundVariables *bg, struct ForegroundVariables *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, FLOAT_P dt_last, bool first_timestep);

FLOAT_P rk2(struct BackgroundVariables *bg, struct ForegroundVariables *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, FLOAT_P dt_last, bool first_timestep);

FLOAT_P rk3(struct BackgroundVariables *bg, struct ForegroundVariables *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, FLOAT_P dt_last, bool first_timestep);

FLOAT_P one_time_step(struct BackgroundVariables *bg, struct ForegroundVariables *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, FLOAT_P dt_last, bool first_timestep);

FLOAT_P get_dt(struct ForegroundVariables *fg, struct GridInfo *grid_info, FLOAT_P dt_last, bool first_timestep);

#endif // ONE_TIME_STEP_H__