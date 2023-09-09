#include "structs.h"
#include "../array_memory_management/array_memory_management.h"
#include <stdlib.h>

void deallocate_foreground_struct_2D(struct ForegroundVariables2D *fg)
{
    deallocate_2D_array(fg->p1);
    deallocate_2D_array(fg->rho1);
    deallocate_2D_array(fg->T1);
    deallocate_2D_array(fg->s1);
    deallocate_2D_array(fg->vx);
    deallocate_2D_array(fg->vz);

    free(fg);
}