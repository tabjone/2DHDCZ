#include <stdlib.h>
#include "background_variables_struct.h"
#include "../../array_utilities/array_memory_management/array_memory_management.h"

void allocate_background_struct(struct BackgroundVariables **bg, int nz_full)
{
    *bg = (struct BackgroundVariables *)malloc(sizeof(struct BackgroundVariables));

    allocate_1D_array(&(*bg)->r, nz_full);
    allocate_1D_array(&(*bg)->T0, nz_full);
    allocate_1D_array(&(*bg)->rho0, nz_full);
    allocate_1D_array(&(*bg)->p0, nz_full);
    allocate_1D_array(&(*bg)->grad_s0, nz_full);
    allocate_1D_array(&(*bg)->g, nz_full);
}