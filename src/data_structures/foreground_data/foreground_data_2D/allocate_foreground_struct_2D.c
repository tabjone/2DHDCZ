#include <stdlib.h>
#include "array_utilities/array_memory_management/array_memory_management.h"
#include "foreground_variables_struct_2D.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"

void allocate_foreground_struct_2D(struct ForegroundVariables2D **fg, struct GridInfo2D *grid_info)
{
    int nz_full = grid_info->nz_full;
    int ny = grid_info->ny;

    
    *fg = (struct ForegroundVariables2D *)malloc(sizeof(struct ForegroundVariables2D));

    allocate_2D_array(&((*fg)->p1), nz_full, ny);
    allocate_2D_array(&((*fg)->T1), nz_full, ny);
    allocate_2D_array(&((*fg)->rho1), nz_full, ny);
    allocate_2D_array(&((*fg)->s1), nz_full, ny);
    allocate_2D_array(&((*fg)->vy), nz_full, ny);
    allocate_2D_array(&((*fg)->vz), nz_full, ny);
}