#ifndef SOLVE_ELLIPTIC_EQUATION_H__
#define SOLVE_ELLIPTIC_EQUATION_H__

#include "shared_files.h"
#include <math.h>

double rhs_elliptic_eq(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, int i, int j);

double rhs_elliptic_eq_vertical_boundary(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg, int i, int j);

void solve_elliptic_equation(struct BackgroundVariables *bg, struct ForegroundVariables2D *fg);

void gauss_seidel(double **A, double *b, double *x, int N, int maxIterations, double tolerance);

#endif // SOLVE_ELLIPTIC_EQUATION_H__