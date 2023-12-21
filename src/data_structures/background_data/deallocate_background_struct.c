#include <stdlib.h>
#include "background_variables_struct.h"
#include "../../array_utilities/array_memory_management/array_memory_management.h"

void deallocate_background_struct(struct BackgroundVariables *bg)
{
    deallocate_1D_array(bg->r);
    deallocate_1D_array(bg->p0);
    deallocate_1D_array(bg->T0);
    deallocate_1D_array(bg->rho0);
    deallocate_1D_array(bg->grad_s0);
    deallocate_1D_array(bg->g);

    free(bg);
}