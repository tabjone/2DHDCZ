#include "structs.h"
#include "../array_memory_management/array_memory_management.h"
#include "global_parameters.h"

void allocate_background_struct(int nz, struct BackgroundVariables *background_variables)
{
    int nz_ghost = UPWIND_ORDER;

    allocate_1D_array(&background_variables->r, nz+2*nz_ghost);
    allocate_1D_array(&background_variables->T0, nz+2*nz_ghost);
    allocate_1D_array(&background_variables->rho0, nz+2*nz_ghost);
    allocate_1D_array(&background_variables->p0, nz+2*nz_ghost);
    allocate_1D_array(&background_variables->grad_s0, nz+2*nz_ghost);
    allocate_1D_array(&background_variables->g, nz+2*nz_ghost);
    background_variables->nz = nz;
    background_variables->nz_ghost = nz_ghost;
    background_variables->nz_full = nz+2*nz_ghost;
}