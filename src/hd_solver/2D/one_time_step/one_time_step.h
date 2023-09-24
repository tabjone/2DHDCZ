#ifndef ONE_TIME_STEP_H__
#define ONE_TIME_STEP_H__

#include "../functions.h"
#include "shared_files.h"

double rk1(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, bool first_run);

double rk2(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, bool first_run);

double one_time_step(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, bool first_run);

#endif // ONE_TIME_STEP_H__