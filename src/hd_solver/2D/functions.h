#ifndef FUNCTIONS_H__
#define FUNCTIONS_H__

#include "./rhs_functions/rhs_functions_2D_hd.h"
#include "./boundaries/boundary_functions_2D_hd.h"
#include "../../shared_files/structs/structs.h"

void one_time_step(struct BackgroundVariables *background_variables, struct ForegroundVariables2D *foreground_variables, int nz, int nx);

#endif // FUNCTIONS_H__