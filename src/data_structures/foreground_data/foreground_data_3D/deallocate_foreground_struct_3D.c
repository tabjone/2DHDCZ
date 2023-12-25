#include <stdlib.h>
#include "array_utilities/array_memory_management/array_memory_management.h"
#include "foreground_variables_struct_3D.h"
#include "data_structures/grid_info/grid_info_3D/grid_info_struct_3D.h"

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