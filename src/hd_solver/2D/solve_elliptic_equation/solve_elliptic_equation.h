#ifndef SOLVE_ELLIPTIC_EQUATION_H__
#define SOLVE_ELLIPTIC_EQUATION_H__

#include "shared_files.h"
#include <math.h>
#include "global_parameters.h"

double rhs_elliptic_eq(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, int i, int j);

double rhs_elliptic_eq_horizontal_boundary(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, struct GridInfo *grid_info, int i, int j);

void solve_elliptic_equation(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg_prev, struct ForegroundVariables2D *fg, struct GridInfo *grid_info);

void gauss_seidel(double **A, double *b, double *x, int N, int maxIterations, double tolerance);

void gauss_seidel_fast_second_order(double **b, double **p1, double **initial_p1, struct GridInfo *grid_info);

#endif // SOLVE_ELLIPTIC_EQUATION_H__