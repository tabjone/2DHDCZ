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

struct ForegroundVariables
{
    FLOAT_P **p1;
    FLOAT_P **rho1;
    FLOAT_P **T1;
    FLOAT_P **s1;
    FLOAT_P **vy;
    FLOAT_P **vz;
};

struct GridInfo
{
    FLOAT_P z_offset;
    int nz, nz_ghost, nz_full;
    FLOAT_P dz, z0, z1;

    int ny;
    FLOAT_P dy, y0, y1;
};

struct MpiInfo
{
    int rank, size;
    bool has_neighbor_below, has_neighbor_above;
    FLOAT_P my_z_offset;
};

void allocate_background_struct(struct BackgroundVariables **bg, struct GridInfo *grid_info);
void allocate_foreground_struct(struct ForegroundVariables **fg, struct GridInfo *grid_info);

void allocate_grid_info_struct(struct GridInfo **grid_info, int nz, int nz_ghost, int nz_full, int ny, FLOAT_P dz, FLOAT_P dy, FLOAT_P z0, FLOAT_P z1, FLOAT_P z_offset, FLOAT_P y0, FLOAT_P y1);

void deallocate_background_struct(struct BackgroundVariables *bg);
void deallocate_foreground_struct(struct ForegroundVariables *fg);
void deallocate_grid_info_struct(struct GridInfo *grid_info);

#endif // STRUCTS_H__