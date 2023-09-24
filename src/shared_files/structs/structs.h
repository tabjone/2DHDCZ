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
};

struct ForegroundVariables2D
{
    double **p1;
    double **rho1;
    double **T1;
    double **s1;
    double **vx;
    double **vz;
};

struct GridInfo
{
    int nz, nz_ghost, nz_full;
    int nx;
    double dz, dx;
    double z0, z1;
    double x0, x1;
};

void allocate_background_struct(struct BackgroundVariables **bg, int nz_full);
void allocate_foreground_struct_2D(struct ForegroundVariables2D **fg, int nz_full, int nx);
void allocate_grid_info_struct(struct GridInfo **grid_info, int nz, int nz_ghost, int nz_full, int nx, double dz, double dx, double z0, double z1, double x0, double x1);

void deallocate_background_struct(struct BackgroundVariables *bg);
void deallocate_foreground_struct_2D(struct ForegroundVariables2D *fg);
void deallocate_grid_info_struct(struct GridInfo *grid_info);

void deep_copy_foreground_2D(struct ForegroundVariables2D *destination, struct ForegroundVariables2D *source, struct GridInfo *grid_info);

#endif // STRUCTS_H__