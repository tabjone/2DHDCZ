#include "structs.h"
#include <stdlib.h>

#if DIMENSIONS == 1
    void allocate_grid_info_struct(struct GridInfo **grid_info, int nz, int nz_ghost, int nz_full, FLOAT_P dz, FLOAT_P z0, FLOAT_P z1)
    {
        // Allocate grid_info
        *grid_info = (struct GridInfo *)malloc(sizeof(struct GridInfo));

        (*grid_info)->z0 = z0;
        (*grid_info)->z1 = z1;
        (*grid_info)->nz = nz;
        (*grid_info)->nz_ghost = nz_ghost;
        (*grid_info)->nz_full = nz_full;
        (*grid_info)->dz = dz;
    }
#elif DIMENSIONS == 2
    void allocate_grid_info_struct(struct GridInfo **grid_info, int nz, int nz_ghost, int nz_full, int ny, FLOAT_P dz, FLOAT_P dy, FLOAT_P z0, FLOAT_P z1, FLOAT_P y0, FLOAT_P y1)
    {
        // Allocate grid_info
        *grid_info = (struct GridInfo *)malloc(sizeof(struct GridInfo));

        (*grid_info)->z0 = z0;
        (*grid_info)->z1 = z1;
        (*grid_info)->y0 = y0;
        (*grid_info)->y1 = y1;
        (*grid_info)->nz = nz;
        (*grid_info)->ny = ny;
        (*grid_info)->nz_ghost = nz_ghost;
        (*grid_info)->nz_full = nz_full;
        (*grid_info)->dy = dy;
        (*grid_info)->dz = dz;
    }
#elif DIMENSIONS == 3
    void allocate_grid_info_struct(struct GridInfo **grid_info, int nz, int nz_ghost, int nz_full, int ny, int nx, FLOAT_P dz, FLOAT_P dy, FLOAT_P dx, FLOAT_P z0, FLOAT_P z1, FLOAT_P y0, FLOAT_P y1, FLOAT_P x0, FLOAT_P x1)
    {
        // Allocate grid_info
        *grid_info = (struct GridInfo *)malloc(sizeof(struct GridInfo));

        (*grid_info)->z0 = z0;
        (*grid_info)->z1 = z1;
        (*grid_info)->x0 = x0;
        (*grid_info)->x1 = x1;
        (*grid_info)->y0 = y0;
        (*grid_info)->y1 = y1;
        (*grid_info)->nz = nz;
        (*grid_info)->nx = nx;
        (*grid_info)->ny = ny;
        (*grid_info)->nz_ghost = nz_ghost;
        (*grid_info)->nz_full = nz_full;
        (*grid_info)->dx = dx;
        (*grid_info)->dz = dz;
        (*grid_info)->dy = dy;
    }
#endif // DIMENSIONS