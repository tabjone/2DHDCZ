#ifndef FUNCTIONS_H__
#define FUNCTIONS_H__

#include "./rhs_functions/rhs_functions.h"
#include "./boundaries/boundaries.h"
#include "./initialization/initialization.h"
#include "shared_files.h"
#include "global_parameters.h"
#include "global_constants.h"

void one_time_step(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg);

#endif // FUNCTIONS_H__