#ifndef EQUATIONS_3D_H__
#define EQUATIONS_3D_H__

#include "shared_files.h"
#include "global_parameters.h"

void equation_of_state_3D(struct ForegroundVariables3D *fg, struct BackgroundVariables *bg, struct GridInfo3D *grid_info);

void first_law_thermodynamics_3D(struct ForegroundVariables3D *fg, struct BackgroundVariables *bg, struct GridInfo3D *grid_info);

#endif // EQUATIONS_3D_H__