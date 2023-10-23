#include "structs.h"
#include <stdlib.h>

void allocate_grid_info_struct_3D(struct GridInfo3D **grid_info, int nz, int nz_ghost, int nz_full, int ny, int nx, FLOAT_P dz, FLOAT_P dy, FLOAT_P dx, FLOAT_P z0, FLOAT_P z1, FLOAT_P z_offset, FLOAT_P y0, FLOAT_P y1, FLOAT_P x0, FLOAT_P x1)
{
    // Allocate grid_info
    *grid_info = (struct GridInfo3D *)malloc(sizeof(struct GridInfo3D));

    (*grid_info)->z0 = z0;
    (*grid_info)->z1 = z1;
    (*grid_info)->y0 = y0;
    (*grid_info)->y1 = y1;
    (*grid_info)->x0 = x0;
    (*grid_info)->x1 = x1;
    (*grid_info)->nx = nx;
    (*grid_info)->nz = nz;
    (*grid_info)->ny = ny;
    (*grid_info)->nz_ghost = nz_ghost;
    (*grid_info)->nz_full = nz_full;
    (*grid_info)->dx = dx;
    (*grid_info)->dy = dy;
    (*grid_info)->dz = dz;
    (*grid_info)->z_offset = z_offset;
}