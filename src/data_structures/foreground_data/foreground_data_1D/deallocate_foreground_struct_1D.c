#include <stdlib.h>
#include "array_utilities/array_memory_management/array_memory_management.h"
#include "foreground_variables_struct_1D.h"
#include "data_structures/grid_info/grid_info_1D/grid_info_struct_1D.h"

void deallocate_foreground_struct_1D(struct ForegroundVariables1D *fg)
{
    deallocate_1D_array(fg->p1);
    deallocate_1D_array(fg->rho1);
    deallocate_1D_array(fg->T1);
    deallocate_1D_array(fg->s1);
    deallocate_1D_array(fg->vz);

    free(fg);
}