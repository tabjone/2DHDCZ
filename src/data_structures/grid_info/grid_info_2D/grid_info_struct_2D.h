#ifndef GRID_INFO_STRUCT_2D_H
#define GRID_INFO_STRUCT_2D_H

#include "global_float_precision.h"

struct GridInfo2D
{
    FLOAT_P z_offset;
    int nz, nz_ghost, nz_full;
    FLOAT_P dz, z0, z1;

    int ny;
    FLOAT_P dy, y0, y1;
};

#endif // GRID_INFO_STRUCT_2D_H