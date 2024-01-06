#ifndef EXTRAPOLATION_2D_H
#define EXTRAPOLATION_2D_H

#include "global_float_precision.h"

void extrapolate_2D_array_antisymmetric_down(FLOAT_P **array, int nz_ghost, int ny);
void extrapolate_2D_array_antisymmetric_up(FLOAT_P **array, int nz_full, int nz_ghost, int ny);
void extrapolate_2D_array_constant_down(FLOAT_P **array, int nz_ghost, int ny);
void extrapolate_2D_array_constant_up(FLOAT_P **array, int nz_full, int nz_ghost, int ny);
void extrapolate_2D_array_symmetric_down(FLOAT_P **array, int nz_ghost, int ny);
void extrapolate_2D_array_symmetric_up(FLOAT_P **array, int nz_full, int nz_ghost, int ny);

#endif // EXTRAPOLATION_2D_H