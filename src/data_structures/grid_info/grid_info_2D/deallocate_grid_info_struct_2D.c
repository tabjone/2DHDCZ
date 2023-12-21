#include <stdlib.h>
#include "grid_info_struct_2D.h"

void deallocate_grid_info_struct_2D(struct GridInfo2D *grid_info)
{
    free(grid_info);
}