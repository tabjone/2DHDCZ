#include "one_time_step.h"

FLOAT_P one_time_step(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, FLOAT_P dt_last, bool first_timestep)
{   
    FLOAT_P dt;
    #if TIME_ORDER == 1
        dt = rk1(bg, fg_prev, fg, grid_info, dt_last, first_timestep);

    #elif TIME_ORDER == 2
        dt = rk2(bg, fg_prev, fg, grid_info, dt_last, first_timestep);
    #elif TIME_ORDER == 3
        dt = rk3(bg, fg_prev, fg, grid_info, dt_last, first_timestep);
    
    #endif
    return dt;
}