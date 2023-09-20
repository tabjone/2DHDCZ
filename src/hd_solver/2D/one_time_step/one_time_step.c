#include "one_time_step.h"

void one_time_step(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, double dt)
{   
    #if TIME_ORDER == 1
        rk1(bg, fg_prev, fg, dt);

    #elif TIME_ORDER == 2
        rk2(bg, fg_prev, fg, dt);
    #endif
}