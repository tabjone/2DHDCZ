#include "structs.h"
#include "../../../shared_files/array_memory_management/array_memory_management.h"

void deallocate_foreground_struct(struct ForegroundVariables *foreground_variables)
{
    deallocate_2D_array(foreground_variables->p1);
    deallocate_2D_array(foreground_variables->rho1);
    deallocate_2D_array(foreground_variables->T1);
    deallocate_2D_array(foreground_variables->s1);
    deallocate_2D_array(foreground_variables->vx);
    deallocate_2D_array(foreground_variables->vz);
}