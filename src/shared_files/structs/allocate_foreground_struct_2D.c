#include "structs.h"
#include "../array_memory_management/array_memory_management.h"
#include "../../global_parameters.h"

void allocate_foreground_struct_2D(int nz, int nx, struct ForegroundVariables2D *foreground_variables)
{
    int nz_ghost = UPWIND_ORDER;

    allocate_2D_array(&foreground_variables->p1, nz+2*nz_ghost, nx);
    allocate_2D_array(&foreground_variables->T1, nz+2*nz_ghost, nx);
    allocate_2D_array(&foreground_variables->rho1, nz+2*nz_ghost, nx);
    allocate_2D_array(&foreground_variables->s1, nz+2*nz_ghost, nx);
    allocate_2D_array(&foreground_variables->vx, nz+2*nz_ghost, nx);
    allocate_2D_array(&foreground_variables->vz, nz+2*nz_ghost, nx);
    foreground_variables->nz = nz;
    foreground_variables->nx = nx;

    foreground_variables->nz_ghost = nz_ghost; 
    foreground_variables->nz_full = nz+2*nz_ghost;
    
}