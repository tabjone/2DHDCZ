#include "one_time_step.h"

FLOAT_P one_time_step(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, bool first_run)
{   
    FLOAT_P dt;
    #if TIME_ORDER == 1
        dt = rk1(bg, fg_prev, fg, grid_info, first_run);

    #elif TIME_ORDER == 2
        dt = rk2(bg, fg_prev, fg, grid_info, first_run);
    #endif
    return dt;
}