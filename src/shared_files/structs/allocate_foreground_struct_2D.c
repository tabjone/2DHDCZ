#include "structs.h"
#include "../array_memory_management/array_memory_management.h"
#include <stdlib.h>

void allocate_foreground_struct_2D(struct ForegroundVariables2D **fg, int nz_full, int nx)
{
    *fg = (struct ForegroundVariables2D *)malloc(sizeof(struct ForegroundVariables2D));

    allocate_2D_array(&((*fg)->p1), nz_full, nx);
    allocate_2D_array(&((*fg)->T1), nz_full, nx);
    allocate_2D_array(&((*fg)->rho1), nz_full, nx);
    allocate_2D_array(&((*fg)->s1), nz_full, nx);
    allocate_2D_array(&((*fg)->vx), nz_full, nx);
    allocate_2D_array(&((*fg)->vz), nz_full, nx);
}