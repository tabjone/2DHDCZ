#include "structs.h"
#include "../../../shared_files/array_memory_management/array_memory_management.h"

void deallocate_background_struct(struct BackgroundVariables *background_variables)
{
    deallocate_1D_array(background_variables->r);
    deallocate_1D_array(background_variables->c_s);
    deallocate_1D_array(background_variables->p0);
    deallocate_1D_array(background_variables->T0);
    deallocate_1D_array(background_variables->rho0);
    deallocate_1D_array(background_variables->Gamma_1);
    deallocate_1D_array(background_variables->H);
    deallocate_1D_array(background_variables->superad_param);
    deallocate_1D_array(background_variables->grad_s0);
    deallocate_1D_array(background_variables->g);
}