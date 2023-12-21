#ifndef GRID_INFO_STRUCT_1D_H
#define GRID_INFO_STRUCT_1D_H

#include "global_float_precision.h"

struct GridInfo1D
{
    FLOAT_P z_offset;
    int nz, nz_ghost, nz_full;
    FLOAT_P dz, z0, z1;
};

#endif // GRID_INFO_STRUCT_1D_H