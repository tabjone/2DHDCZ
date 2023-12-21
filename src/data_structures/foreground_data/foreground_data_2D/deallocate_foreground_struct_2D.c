#include <stdlib.h>
#include "../../array_utilities/array_memory_management/array_memory_management.h"
#include "foreground_variables_struct_2D.h"
#include "../grid_info/grid_info_2D/grid_info_struct_2D.h"

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