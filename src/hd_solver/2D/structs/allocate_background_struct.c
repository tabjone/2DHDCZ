#include "structs.h"
#include "../../../shared_files/array_memory_management/array_memory_management.h"

void allocate_background_struct(int nz, struct BackgroundVariables *background_variables)
{
    allocate_1D_array(&background_variables->r, nz);
    allocate_1D_array(&background_variables->T0, nz);
    allocate_1D_array(&background_variables->rho0, nz);
    allocate_1D_array(&background_variables->p0, nz);
    allocate_1D_array(&background_variables->c_s, nz);
    allocate_1D_array(&background_variables->Gamma_1, nz);
    allocate_1D_array(&background_variables->H, nz);
    allocate_1D_array(&background_variables->superad_param, nz);
    allocate_1D_array(&background_variables->grad_s0, nz);
    allocate_1D_array(&background_variables->g, nz);
}