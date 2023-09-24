#include "structs.h"
#include <stdlib.h>

void deallocate_grid_info_struct(struct GridInfo *grid_info)
{
    free(grid_info);
}