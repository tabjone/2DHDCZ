#ifndef ONE_TIME_STEP_H__
#define ONE_TIME_STEP_H__

#include "../functions.h"
#include "shared_files.h"

void rk1(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, double dt);


void one_time_step(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, double dt);

#endif // ONE_TIME_STEP_H__