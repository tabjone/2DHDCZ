#include "structs.h"
#include "../array_memory_management/array_memory_management.h"
#include "global_parameters.h"
#include <stdlib.h>

void allocate_background_struct(int nz, struct BackgroundVariables **bg)
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

    *bg = (struct BackgroundVariables *)malloc(sizeof(struct BackgroundVariables));

    allocate_1D_array(&(*bg)->r, nz+2*nz_ghost);
    allocate_1D_array(&(*bg)->T0, nz+2*nz_ghost);
    allocate_1D_array(&(*bg)->rho0, nz+2*nz_ghost);
    allocate_1D_array(&(*bg)->p0, nz+2*nz_ghost);
    allocate_1D_array(&(*bg)->grad_s0, nz+2*nz_ghost);
    allocate_1D_array(&(*bg)->g, nz+2*nz_ghost);
    (*bg)->nz = nz;
    (*bg)->nz_ghost = nz_ghost;
    (*bg)->nz_full = nz+2*nz_ghost;
}