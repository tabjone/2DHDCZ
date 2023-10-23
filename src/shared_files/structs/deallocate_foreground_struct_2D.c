#include "structs.h"
#include <stdlib.h>

void deallocate_foreground_struct_2D(struct ForegroundVariables2D *fg)
{
    deallocate_2D_array(fg->p1);
    deallocate_2D_array(fg->rho1);
    deallocate_2D_array(fg->T1);
    deallocate_2D_array(fg->s1);
    deallocate_2D_array(fg->vz);
    deallocate_2D_array(fg->vy);

    free(fg);
}