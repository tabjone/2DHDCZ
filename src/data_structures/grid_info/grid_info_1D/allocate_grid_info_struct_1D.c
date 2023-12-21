#include <stdlib.h>
#include "grid_info_struct_1D.h"

void allocate_grid_info_struct_1D(struct GridInfo1D **grid_info)
{
    // Allocate grid_info
    *grid_info = (struct GridInfo1D *)malloc(sizeof(struct GridInfo1D));
}