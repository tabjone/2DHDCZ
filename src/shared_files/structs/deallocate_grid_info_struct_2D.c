#include "structs.h"
#include <stdlib.h>

void deallocate_grid_info_struct_2D(struct GridInfo2D *grid_info)
{
    free(grid_info);
}