#ifndef DERIVATIVES_3D_H
#define DERIVATIVES_3D_H

#include "global_float_precision.h"

FLOAT_P central_first_derivative_x_3D(FLOAT_P ***array, int i, int j, int k, int nx, FLOAT_P one_over_2dx);
FLOAT_P central_first_derivative_y_3D(FLOAT_P ***array, int i, int j, int k, int ny, FLOAT_P one_over_2dy);
FLOAT_P central_first_derivative_z_3D(FLOAT_P ***array, int i, int j, int k, FLOAT_P one_over_2dz);

FLOAT_P central_second_derivative_x_3D(FLOAT_P ***array, int i, int j, int k, int nx, FLOAT_P one_over_dxdx);
FLOAT_P central_second_derivative_y_3D(FLOAT_P ***array, int i, int j, int k, int ny, FLOAT_P one_over_dydy);
FLOAT_P central_second_derivative_z_3D(FLOAT_P ***array, int i, int j, int k, FLOAT_P one_over_dzdz);

FLOAT_P upwind_first_derivative_x_3D(FLOAT_P ***array, FLOAT_P ***velocity, int i, int j, int k, int nx, FLOAT_P one_over_dx, FLOAT_P one_over_2dx);
FLOAT_P upwind_first_derivative_y_3D(FLOAT_P ***array, FLOAT_P ***velocity, int i, int j, int k, int ny, FLOAT_P one_over_dy, FLOAT_P one_over_2dy);
FLOAT_P upwind_first_derivative_z_3D(FLOAT_P ***array, FLOAT_P ***velocity, int i, int j, int k, int nz, FLOAT_P one_over_dz, FLOAT_P one_over_2dz);

#endif // DERIVATIVES_3D_H