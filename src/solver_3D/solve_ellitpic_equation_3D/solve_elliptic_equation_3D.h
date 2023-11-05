#ifndef SOLVE_ELLIPTIC_EQUATION_H__
#define SOLVE_ELLIPTIC_EQUATION_H__

#include "shared_files.h"
#include "global_parameters.h"
#include "../rhs_functions_3D/rhs_functions_3D.h"

void solve_elliptic_equation_3D(struct BackgroundVariables *bg, struct ForegroundVariables3D *fg_prev, struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info);

void gauss_seidel_3D(FLOAT_P ***b, FLOAT_P ***p1, FLOAT_P ***initial_p1, struct GridInfo3D *grid_info);
void gauss_seidel_3D_mpi(FLOAT_P ***b, FLOAT_P ***p1, FLOAT_P ***initial_p1, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info);

void communicate_p_gauss_seidel_3D(FLOAT_P ***array, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info);

#endif // SOLVE_ELLIPTIC_EQUATION_H__