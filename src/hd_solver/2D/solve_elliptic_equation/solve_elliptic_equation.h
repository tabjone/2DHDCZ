#ifndef SOLVE_ELLIPTIC_EQUATION_H__
#define SOLVE_ELLIPTIC_EQUATION_H__

#include "shared_files.h"
#include <math.h>

double rhs_elliptic_eq(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, int i, int j);

double rhs_elliptic_eq_horizontal_boundary(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, int i, int j);

void solve_elliptic_equation(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg);

void gauss_seidel(double **A, double *b, double *x, int N, int maxIterations, double tolerance);

void gauss_seidel_fast_second_order(double **b, double dz, double dx, int nx, int nz, int nz_ghost, int max_iterations, double tolerance, double **p1);

#endif // SOLVE_ELLIPTIC_EQUATION_H__