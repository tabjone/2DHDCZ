#ifndef EQUATIONS_H__
#define EQUATIONS_H__

#include "shared_files.h"
#include "global_parameters.h"
#include "global_parameters.h"

void equation_of_state(struct ForegroundVariables *fg, struct BackgroundVariables *bg, struct GridInfo *grid_info);
void first_law_thermodynamics(struct ForegroundVariables *fg, struct BackgroundVariables *bg, struct GridInfo *grid_info);

#endif // EQUATIONS_H__