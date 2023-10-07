#include "structs.h"
#include "../array_memory_management/array_memory_management.h"
#include <stdlib.h>

void deallocate_background_struct(struct BackgroundVariables *bg)
{
    deallocate_1D_array(bg->r);
    deallocate_1D_array(bg->p0);
    deallocate_1D_array(bg->T0);
    deallocate_1D_array(bg->rho0);
    deallocate_1D_array(bg->one_over_rho0);
    deallocate_1D_array(bg->grad_s0);
    deallocate_1D_array(bg->g);
    deallocate_1D_array(bg->grad_g);
    deallocate_1D_array(bg->grad_rho0);

    free(bg);
}