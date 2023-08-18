#ifndef hd_2D_functions_H__
#define hd_2D_functions_H__

#include "./rhs_functions/rhs_functions_2D_hd.h"
#include "./boundaries/boundary_functions_2D_hd.h"

int main_hd_2D(int argc, char *argv[]);
void one_time_step_2D_hd(double **s1, double **p1, double **rho1, double **T1, double **vx, double **vz, double *grad_s0, double *p0, double *rho0, double *T0, double *g, int nz, int nx, double dx, double dz);

#endif // hd_2D_functions_H__