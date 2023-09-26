#ifndef ONE_TIME_STEP_H__
#define ONE_TIME_STEP_H__

#include "../functions.h"
#include "shared_files.h"

FLOAT_P rk1(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, bool first_run);

FLOAT_P rk2(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, bool first_run);

FLOAT_P one_time_step(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, bool first_run);

#endif // ONE_TIME_STEP_H__