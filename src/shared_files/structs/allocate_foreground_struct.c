#include "structs.h"
#include "../array_memory_management/array_memory_management.h"
#include <stdlib.h>

void allocate_foreground_struct(struct ForegroundVariables **fg, struct GridInfo *grid_info)
{
    int nz_full = grid_info->nz_full;
    
    *fg = (struct ForegroundVariables *)malloc(sizeof(struct ForegroundVariables));

    #if DIMENSIONS == 1
        allocate_1D_array(&((*fg)->p1), nz_full);
        allocate_1D_array(&((*fg)->T1), nz_full);
        allocate_1D_array(&((*fg)->rho1), nz_full);
        allocate_1D_array(&((*fg)->s1), nz_full);
        allocate_1D_array(&((*fg)->vz), nz_full);
        #if BFIELD_ON == 1
            allocate_1D_array(&((*fg)->Bz), nz_full);
        #endif // BFIELD_ON
    #elif DIMENSIONS == 2
        int ny = grid_info->ny;
        allocate_2D_array(&((*fg)->p1), nz_full, ny);
        allocate_2D_array(&((*fg)->T1), nz_full, ny);
        allocate_2D_array(&((*fg)->rho1), nz_full, ny);
        allocate_2D_array(&((*fg)->s1), nz_full, ny);
        allocate_2D_array(&((*fg)->vy), nz_full, ny);
        allocate_2D_array(&((*fg)->vz), nz_full, ny);
        #if BFIELD_ON == 1
            allocate_2D_array(&((*fg)->Bz), nz_full, ny);
            allocate_2D_array(&((*fg)->By), nz_full, ny);
        #endif // BFIELD_ON
    #elif DIMENSIONS == 3
        int nx = grid_info->nx;
        int ny = grid_info->ny;
        allocate_3D_array(&((*fg)->p1), nz_full, ny, nx);
        allocate_3D_array(&((*fg)->T1), nz_full, ny, nx);
        allocate_3D_array(&((*fg)->rho1), nz_full, ny, nx);
        allocate_3D_array(&((*fg)->s1), nz_full, ny, nx);
        allocate_3D_array(&((*fg)->vx), nz_full, ny, nx);
        allocate_3D_array(&((*fg)->vy), nz_full, ny, nx);
        allocate_3D_array(&((*fg)->vz), nz_full, ny, nx);
        #if BFIELD_ON == 1
            allocate_3D_array(&((*fg)->Bz), nz_full, ny, nx);
            allocate_3D_array(&((*fg)->Bx), nz_full, ny, nx);
            allocate_3D_array(&((*fg)->By), nz_full, ny, nx);
        #endif // BFIELD_ON
    #endif // DIMENSIONS
}