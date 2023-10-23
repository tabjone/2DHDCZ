#ifndef SOLVE_ELLIPTIC_EQUATION_H__
#define SOLVE_ELLIPTIC_EQUATION_H__

#include "shared_files.h"
#include <math.h>
#include "global_parameters.h"
#include "../rhs_functions/rhs_functions.h"

void gauss_seidel_2D(FLOAT_P **b, FLOAT_P **p1, FLOAT_P **initial_p1, struct GridInfo2D *grid_info);

void solve_elliptic_equation(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info);

void communicate_p_gauss_seidel(FLOAT_P **array, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info);
void gauss_seidel_2D_mpi(FLOAT_P **b, FLOAT_P **p1, FLOAT_P **initial_p1, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info);

#endif // SOLVE_ELLIPTIC_EQUATION_H__