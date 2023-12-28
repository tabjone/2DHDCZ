#ifndef DERIVATIVES_1D_H
#define DERIVATIVES_1D_H

#include "global_float_precision.h"

FLOAT_P central_first_derivative_z_1D(FLOAT_P *array, int i, FLOAT_P one_over_2dz);
FLOAT_P central_second_derivative_z_1D(FLOAT_P *array, int i, FLOAT_P one_over_dzdz);@
FLOAT_P upwind_first_derivative_z_1D(FLOAT_P *array, FLOAT_P *velocity, int i, FLOAT_P one_over_dz, FLOAT_P one_over_2dz);

#endif // DERIVATIVES_1D_H