#include "structs.h"
#include <stdlib.h>

void deallocate_grid_info_struct_3D(struct GridInfo3D *grid_info)
{
    free(grid_info);
}