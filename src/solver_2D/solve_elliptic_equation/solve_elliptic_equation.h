#ifndef SOLVE_ELLIPTIC_EQUATION_H__
#define SOLVE_ELLIPTIC_EQUATION_H__

#include "shared_files.h"
#include <math.h>
#include "global_parameters.h"
#include "../rhs_functions/rhs_functions.h"
#include "global_boundary.h"
#include "../../iterative_solver_module/iterative_solver_module.h"


void solve_elliptic_equation(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables *precalc);


#endif // SOLVE_ELLIPTIC_EQUATION_H__