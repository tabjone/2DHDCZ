#ifndef SPACIAL_DERIVATIVES_3D_H__
#define SPACIAL_DERIVATIVES_3D_H__

#include "shared_files.h"
#include "global_parameters.h"

FLOAT_P central_first_derivative_x_3D(FLOAT_P ***array, int i, int j, int k, FLOAT_P dx, int nx);
FLOAT_P central_first_derivative_y_3D(FLOAT_P ***array, int i, int j, int k, FLOAT_P dy, int ny);
FLOAT_P central_first_derivative_z_3D(FLOAT_P ***array, int i, int j, int k, FLOAT_P dz, int nz);

FLOAT_P central_second_derivative_x_3D(FLOAT_P ***array, int i, int j, int k, FLOAT_P dx, int nx);
FLOAT_P central_second_derivative_y_3D(FLOAT_P ***array, int i, int j, int k, FLOAT_P dy, int ny);
FLOAT_P central_second_derivative_z_3D(FLOAT_P ***array, int i, int j, int k, FLOAT_P dz, int nz);

FLOAT_P upwind_first_derivative_x_3D(FLOAT_P ***array, FLOAT_P ***velocity, int i, int j, int k, FLOAT_P dx, int nx);
FLOAT_P upwind_first_derivative_y_3D(FLOAT_P ***array, FLOAT_P ***velocity, int i, int j, int k, FLOAT_P dy, int ny);
FLOAT_P upwind_first_derivative_z_3D(FLOAT_P ***array, FLOAT_P ***velocity, int i, int j, int k, FLOAT_P dz, int nz);

#endif // SPACIAL_DERIVATIVES_3D_H__