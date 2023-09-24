#include "structs.h"
#include <stdlib.h>

void allocate_grid_info_struct(struct GridInfo **grid_info, int nz, int nz_ghost, int nz_full, int nx, double dz, double dx, double z0, double z1, double x0, double x1)
{
    // Allocate grid_info
    *grid_info = (struct GridInfo *)malloc(sizeof(struct GridInfo));

    (*grid_info)->z0 = z0;
    (*grid_info)->z1 = z1;
    (*grid_info)->x0 = x0;
    (*grid_info)->x1 = x1;
    (*grid_info)->nz = nz;
    (*grid_info)->nx = nx;
    (*grid_info)->nz_ghost = nz_ghost;
    (*grid_info)->nz_full = nz_full;
    (*grid_info)->dx = dx;
    (*grid_info)->dz = dz;

}