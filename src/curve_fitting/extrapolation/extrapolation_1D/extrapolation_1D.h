#ifndef EXTRAPOLATION_1D_H
#define EXTRAPOLATION_1D_H

#include "global_float_precision.h"

void extrapolate_1D_array_constant_down(FLOAT_P *array, int nz_ghost);
void extrapolate_1D_array_constant_up(FLOAT_P *array, int nz_full, int nz_ghost);

#endif // EXTRAPOLATION_1D_H