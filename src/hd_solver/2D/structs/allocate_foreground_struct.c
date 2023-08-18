#include "structs.h"
#include "../../../shared_files/array_memory_management/array_memory_management.h"

void allocate_foreground_struct(int nz, int nx, struct ForegroundVariables *foreground_variables)
{
    allocate_2D_array(&foreground_variables->p1, nz, nx);
    allocate_2D_array(&foreground_variables->T1, nz, nx);
    allocate_2D_array(&foreground_variables->rho1, nz, nx);
    allocate_2D_array(&foreground_variables->s1, nz, nx);
    allocate_2D_array(&foreground_variables->vx, nz, nx);
    allocate_2D_array(&foreground_variables->vz, nz, nx);
}