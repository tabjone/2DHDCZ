#ifndef RHS_FUNCTIONS_2D_HD_H
#define RHS_FUNCTIONS_2D_HD_H

double rhs_dvx_dt_2D_hd(double **p1, double **vx, double **vz, double *rho0, int i, int j, double dx, double dz);

double rhs_dvz_dt_2D_hd(double **rho1, double **p1, double **vx, double **vz, double *rho0, double *g, int i, int j, double dx, double dz);

double rhs_ds1_dt_2D_hd(double **s1, double *s0, double **vx, double **vz, double *T0, double *rho0, int i, int j, double dx, double dz);

#endif // RHS_FUNCTIONS_2D_HD_H