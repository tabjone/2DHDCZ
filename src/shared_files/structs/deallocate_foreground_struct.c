#include "structs.h"
#include "../array_memory_management/array_memory_management.h"
#include <stdlib.h>

void deallocate_foreground_struct(struct ForegroundVariables *fg)
{
    deallocate_2D_array(fg->p1);
    deallocate_2D_array(fg->rho1);
    deallocate_2D_array(fg->T1);
    deallocate_2D_array(fg->s1);
    deallocate_2D_array(fg->vz);
    deallocate_2D_array(fg->vy);

    free(fg);
}