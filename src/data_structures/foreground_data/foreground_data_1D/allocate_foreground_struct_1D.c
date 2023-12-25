#include <stdlib.h>
#include "array_utilities/array_memory_management/array_memory_management.h"
#include "foreground_variables_struct_1D.h"
#include "data_structures/grid_info/grid_info_1D/grid_info_struct_1D.h"

void allocate_foreground_struct_1D(struct ForegroundVariables1D **fg, struct GridInfo1D *grid_info)
{
    int nz_full = grid_info->nz_full;
    int ny = grid_info->ny;

    
    *fg = (struct ForegroundVariables1D *)malloc(sizeof(struct ForegroundVariables1D));

    allocate_1D_array(&((*fg)->p1), nz_full);
    allocate_1D_array(&((*fg)->T1), nz_full);
    allocate_1D_array(&((*fg)->rho1), nz_full);
    allocate_1D_array(&((*fg)->s1), nz_full);
    allocate_1D_array(&((*fg)->vz), nz_full);
}