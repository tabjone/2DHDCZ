#include <stdlib.h>
#include "grid_info_struct_3D.h"

void allocate_grid_info_struct_3D(struct GridInfo3D **grid_info)
{
    // Allocate grid_info
    *grid_info = (struct GridInfo3D *)malloc(sizeof(struct GridInfo3D));
}