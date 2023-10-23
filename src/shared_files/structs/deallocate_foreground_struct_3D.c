#include "structs.h"
#include <stdlib.h>

void deallocate_foreground_struct_3D(struct ForegroundVariables3D *fg)
{
    deallocate_3D_array(fg->p1);
    deallocate_3D_array(fg->rho1);
    deallocate_3D_array(fg->T1);
    deallocate_3D_array(fg->s1);
    deallocate_3D_array(fg->vz);
    deallocate_3D_array(fg->vy);
    deallocate_3D_array(fg->vx);

    free(fg);
}