#include "structs.h"
#include "../array_memory_management/array_memory_management.h"
#include "global_parameters.h"
#include <stdlib.h>

void allocate_foreground_struct_2D(int nz, int nx, double dz, double dx, struct ForegroundVariables2D **fg)
{
    int nz_ghost;
    if (UPWIND_ORDER >= CENTRAL_ORDER)
    {
        nz_ghost = UPWIND_ORDER;
    }
    else
    {
        nz_ghost = CENTRAL_ORDER;
    }

    *fg = (struct ForegroundVariables2D *)malloc(sizeof(struct ForegroundVariables2D));

    allocate_2D_array(&((*fg)->p1), nz+2*nz_ghost, nx);
    allocate_2D_array(&((*fg)->T1), nz+2*nz_ghost, nx);
    allocate_2D_array(&((*fg)->rho1), nz+2*nz_ghost, nx);
    allocate_2D_array(&((*fg)->s1), nz+2*nz_ghost, nx);
    allocate_2D_array(&((*fg)->vx), nz+2*nz_ghost, nx);
    allocate_2D_array(&((*fg)->vz), nz+2*nz_ghost, nx);

    (*fg)->nz = nz;
    (*fg)->nx = nx;
    (*fg)->nz_ghost = nz_ghost;
    (*fg)->nz_full = nz+2*nz_ghost;
    (*fg)->dx = dx;
    (*fg)->dz = dz;
}