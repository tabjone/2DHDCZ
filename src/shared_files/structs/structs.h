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
    FLOAT_P *eta_over_four_pi_rho0_T0;
};

struct ForegroundVariables2D
{
    FLOAT_P **p1;
    FLOAT_P **rho1;
    FLOAT_P **T1;
    FLOAT_P **s1;
    FLOAT_P **vy;
    FLOAT_P **vz;
};

struct ForegroundVariables3D
{
    FLOAT_P ***p1;
    FLOAT_P ***rho1;
    FLOAT_P ***T1;
    FLOAT_P ***s1;
    FLOAT_P ***vx;
    FLOAT_P ***vy;
    FLOAT_P ***vz;
};

struct GridInfo2D
{
    FLOAT_P z_offset;
    int nz, nz_ghost, nz_full;
    FLOAT_P dz, z0, z1;

    int ny;
    FLOAT_P dy, y0, y1;
};

struct GridInfo3D
{
    FLOAT_P z_offset;
    int nz, nz_ghost, nz_full;
    FLOAT_P dz, z0, z1;

    int ny;
    FLOAT_P dy, y0, y1;

    int nx;
    FLOAT_P dx, x0, x1;
};

struct MpiInfo
{
    int rank, size;
    bool has_neighbor_below, has_neighbor_above;
    FLOAT_P my_z_offset;
};

#include "../array_memory_management/array_memory_management.h"

void allocate_background_struct(struct BackgroundVariables **bg, int nz_full);
void deallocate_background_struct(struct BackgroundVariables *bg);

void allocate_foreground_struct_2D(struct ForegroundVariables2D **fg, struct GridInfo2D *grid_info);
void allocate_foreground_struct_3D(struct ForegroundVariables3D **fg, struct GridInfo3D *grid_info);

void deallocate_foreground_struct_2D(struct ForegroundVariables2D *fg);
void deallocate_foreground_struct_3D(struct ForegroundVariables3D *fg);

void allocate_grid_info_struct_2D(struct GridInfo2D **grid_info, int nz, int nz_ghost, int nz_full, int ny, FLOAT_P dz, FLOAT_P dy, FLOAT_P z0, FLOAT_P z1, FLOAT_P z_offset, FLOAT_P y0, FLOAT_P y1);
void allocate_grid_info_struct_3D(struct GridInfo3D **grid_info, int nz, int nz_ghost, int nz_full, int ny, int nx, FLOAT_P dz, FLOAT_P dy, FLOAT_P dx, FLOAT_P z0, FLOAT_P z1, FLOAT_P z_offset, FLOAT_P y0, FLOAT_P y1, FLOAT_P x0, FLOAT_P x1);

void deallocate_grid_info_struct_2D(struct GridInfo2D *grid_info);
void deallocate_grid_info_struct_3D(struct GridInfo3D *grid_info);

#endif // STRUCTS_H__