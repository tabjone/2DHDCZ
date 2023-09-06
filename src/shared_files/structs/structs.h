#ifndef STRUCTS_H__
#define STRUCTS_H__

struct BackgroundVariables
{
    double *r;
    double *p0;
    double *T0;
    double *rho0;
    double *grad_s0;
    double *g;
    int nz, nz_ghost, nz_full;
};

struct ForegroundVariables2D
{
    double **p1;
    double **rho1;
    double **T1;
    double **s1;
    double **vx;
    double **vz;
    int nz, nz_ghost, nz_full;
    int nx;
};

void allocate_background_struct(int nz, struct BackgroundVariables *background_variables);
void allocate_foreground_struct_2D(int nz, int nx, struct ForegroundVariables2D *foreground_variables);

void deallocate_background_struct(struct BackgroundVariables *background_variables);
void deallocate_foreground_struct_2D(struct ForegroundVariables2D *foreground_variables);

#endif // STRUCTS_H__