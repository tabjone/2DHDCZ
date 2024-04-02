#ifndef DERIVATIVES_2D_H
#define DERIVATIVES_2D_H

#include "global_float_precision.h"

FLOAT_P central_first_derivative_y_2D(FLOAT_P **array, int i, int j, int ny, FLOAT_P one_over_2dy);
FLOAT_P central_first_derivative_z_2D(FLOAT_P **array, int i, int j, FLOAT_P one_over_2dz);

FLOAT_P central_second_derivative_y_2D(FLOAT_P **array, int i, int j, int ny, FLOAT_P one_over_dydy);
FLOAT_P central_second_derivative_yz_2D(FLOAT_P **array, int i, int j, int ny, FLOAT_P one_over_4dydz);
FLOAT_P central_second_derivative_z_2D(FLOAT_P **array, int i, int j, FLOAT_P one_over_dzdz);

FLOAT_P central_third_derivative_yyz_2D(FLOAT_P **array, int i, int j, int ny, FLOAT_P one_over_8dydydz);
FLOAT_P central_third_derivative_yzz_2D(FLOAT_P **array, int i, int j, int ny, FLOAT_P one_over_8dydzdz);

FLOAT_P central_third_derivative_z_2D(FLOAT_P **array, int i, int j, FLOAT_P one_over_2dzdzdz);
FLOAT_P central_third_derivative_y_2D(FLOAT_P **array, int i, int j, int ny, FLOAT_P one_over_2dydydy);

FLOAT_P upwind_first_derivative_y_2D(FLOAT_P **array, FLOAT_P **velocity, int i, int j, int ny, FLOAT_P one_over_dy, FLOAT_P one_over_2dy);
FLOAT_P upwind_first_derivative_z_2D(FLOAT_P **array, FLOAT_P **velocity, int i, int j, FLOAT_P one_over_dz, FLOAT_P one_over_2dz);

FLOAT_P upwind_y_second_derivative_yz_2D(FLOAT_P **array, FLOAT_P **velocity, int i, int j, int ny, FLOAT_P one_over_2dydz);
FLOAT_P upwind_z_second_derivative_yz_2D(FLOAT_P **array, FLOAT_P **velocity, int i, int j, int ny, FLOAT_P one_over_2dydz);

#endif // DERIVATIVES_2D_H