#ifndef RHS_FUNCTIONS_H__
#define RHS_FUNCTIONS_H__

#include "shared_files.h"
#include "global_parameters.h"
#include "../boundaries/boundaries.h"

double rhs_dvx_dt(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, int i, int j);

double rhs_dvz_dt(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, int i, int j);

double rhs_ds1_dt(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, int i, int j);

#endif // RHS_FUNCTIONS_H__