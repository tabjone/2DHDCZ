#include "structs.h"
#include "../array_memory_management/array_memory_management.h"
#include <stdlib.h>

void allocate_foreground_struct(struct ForegroundVariables **fg, struct GridInfo *grid_info)
{
    int nz_full = grid_info->nz_full;
    int ny = grid_info->ny;

    
    *fg = (struct ForegroundVariables *)malloc(sizeof(struct ForegroundVariables));

    allocate_2D_array(&((*fg)->p1), nz_full, ny);
    allocate_2D_array(&((*fg)->T1), nz_full, ny);
    allocate_2D_array(&((*fg)->rho1), nz_full, ny);
    allocate_2D_array(&((*fg)->s1), nz_full, ny);
    allocate_2D_array(&((*fg)->vy), nz_full, ny);
    allocate_2D_array(&((*fg)->vz), nz_full, ny);
}