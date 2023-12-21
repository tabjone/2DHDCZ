#include <stdlib.h>
#include "grid_info_struct_1D.h"

void deallocate_grid_info_struct_1D(struct GridInfo1D *grid_info)
{
    free(grid_info);
}