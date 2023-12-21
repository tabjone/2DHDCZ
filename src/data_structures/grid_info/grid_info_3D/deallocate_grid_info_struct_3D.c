#include <stdlib.h>
#include "grid_info_struct_3D.h"

void deallocate_grid_info_struct_3D(struct GridInfo3D *grid_info)
{
    free(grid_info);
}