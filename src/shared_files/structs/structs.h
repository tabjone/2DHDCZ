#ifndef STRUCTS_H__
#define STRUCTS_H__

#include "global_parameters.h"

struct BackgroundVariables
{
    FLOAT_P *r;
    FLOAT_P *p0;
    FLOAT_P *T0;
    FLOAT_P *rho0;
    FLOAT_P *one_over_rho0;
    FLOAT_P *grad_s0;
    FLOAT_P *g;
    FLOAT_P *grad_g;
    FLOAT_P *grad_rho0;
};

struct ForegroundVariables
{
    #if DIMENSIONS == 1
        FLOAT_P *p1;
        FLOAT_P *rho1;
        FLOAT_P *T1;
        FLOAT_P *s1;
        FLOAT_P *vz;
        FLOAT_P *Bz;
    #elif DIMENSIONS == 2
        FLOAT_P **p1;
        FLOAT_P **rho1;
        FLOAT_P **T1;
        FLOAT_P **s1;
        FLOAT_P **vy;
        FLOAT_P **vz;
        FLOAT_P **By;
        FLOAT_P **Bz;
    #elif DIMENSIONS == 3
        FLOAT_P ***p1;
        FLOAT_P ***rho1;
        FLOAT_P ***T1;
        FLOAT_P ***s1;
        FLOAT_P ***vx;
        FLOAT_P ***vy;
        FLOAT_P ***vz;
        FLOAT_P ***Bx;
        FLOAT_P ***By;
        FLOAT_P ***Bz;
    #endif // DIMENSIONS
};

struct GridInfo
{
    int nz, nz_ghost, nz_full;
    FLOAT_P dz, z0, z1;
    #if DIMENSIONS == 2
        int ny;
        FLOAT_P dy, y0, y1;
    #elif DIMENSIONS == 3
        int ny, nx;
        FLOAT_P dy, dx, y0, y1, x0, x1;
    #endif // DIMENSIONS
};

struct MpiInfo
{
    int rank, size;  
};

void allocate_background_struct(struct BackgroundVariables **bg, struct GridInfo *grid_info);
void allocate_foreground_struct(struct ForegroundVariables **fg, struct GridInfo *grid_info);

#if DIMENSIONS == 1
    void allocate_grid_info_struct(struct GridInfo **grid_info, int nz, int nz_ghost, int nz_full, FLOAT_P dz, FLOAT_P z0, FLOAT_P z1);
#elif DIMENSIONS == 2
    void allocate_grid_info_struct(struct GridInfo **grid_info, int nz, int nz_ghost, int nz_full, int nx, FLOAT_P dz, FLOAT_P dy, FLOAT_P z0, FLOAT_P z1, FLOAT_P y0, FLOAT_P y1);
#elif DIMENSIONS == 3
    void allocate_grid_info_struct(struct GridInfo **grid_info, int nz, int nz_ghost, int nz_full, int ny, int nx, FLOAT_P dz, FLOAT_P dy, FLOAT_P dx, FLOAT_P z0, FLOAT_P z1, FLOAT_P y0, FLOAT_P y1, FLOAT_P x0, FLOAT_P x1);
#endif // DIMENSIONS

void deallocate_background_struct(struct BackgroundVariables *bg);
void deallocate_foreground_struct(struct ForegroundVariables *fg);
void deallocate_grid_info_struct(struct GridInfo *grid_info);

#endif // STRUCTS_H__