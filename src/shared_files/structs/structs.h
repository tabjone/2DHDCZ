#ifndef STRUCTS_H__
#define STRUCTS_H__

#include "global_parameters.h"

struct BackgroundVariables
{
    FLOAT_P *r;
    FLOAT_P *p0;
    FLOAT_P *T0;
    FLOAT_P *rho0;
    FLOAT_P *grad_s0;
    FLOAT_P *g;
};

struct ForegroundVariables2D
{
    FLOAT_P **p1;
    FLOAT_P **rho1;
    FLOAT_P **T1;
    FLOAT_P **s1;
    FLOAT_P **vx;
    FLOAT_P **vz;
};

struct GridInfo
{
    int nz, nz_ghost, nz_full;
    int nx;
    FLOAT_P dz, dx;
    FLOAT_P z0, z1;
    FLOAT_P x0, x1;
};

void allocate_background_struct(struct BackgroundVariables **bg, int nz_full);
void allocate_foreground_struct_2D(struct ForegroundVariables2D **fg, int nz_full, int nx);
void allocate_grid_info_struct(struct GridInfo **grid_info, int nz, int nz_ghost, int nz_full, int nx, FLOAT_P dz, FLOAT_P dx, FLOAT_P z0, FLOAT_P z1, FLOAT_P x0, FLOAT_P x1);

void deallocate_background_struct(struct BackgroundVariables *bg);
void deallocate_foreground_struct_2D(struct ForegroundVariables2D *fg);
void deallocate_grid_info_struct(struct GridInfo *grid_info);

void deep_copy_foreground_2D(struct ForegroundVariables2D *destination, struct ForegroundVariables2D *source, struct GridInfo *grid_info);

#endif // STRUCTS_H__