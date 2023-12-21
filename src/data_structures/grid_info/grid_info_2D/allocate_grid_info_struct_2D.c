#include <stdlib.h>
#include "grid_info_struct_2D.h"

void allocate_grid_info_struct_2D(struct GridInfo2D **grid_info)
{
    // Allocate grid_info
    *grid_info = (struct GridInfo2D *)malloc(sizeof(struct GridInfo2D));
}