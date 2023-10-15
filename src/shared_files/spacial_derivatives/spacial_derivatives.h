#ifndef SPACIAL_DERIVATIVES_H__
#define SPACIAL_DERIVATIVES_H__

#include "global_parameters.h"
#include "shared_files.h"

FLOAT_P backward_first_derivative_first_order(FLOAT_P centre, FLOAT_P left, FLOAT_P dx_);
FLOAT_P backward_first_derivative_second_order(FLOAT_P centre, FLOAT_P left, FLOAT_P left2, FLOAT_P dx_);

FLOAT_P central_first_derivative_second_order(FLOAT_P left, FLOAT_P right, FLOAT_P dx_);
FLOAT_P central_second_derivative_second_order(FLOAT_P centre, FLOAT_P left, FLOAT_P right, FLOAT_P dx_);

FLOAT_P forward_first_derivative_first_order(FLOAT_P centre, FLOAT_P right, FLOAT_P dx_);
FLOAT_P forward_first_derivative_second_order(FLOAT_P centre, FLOAT_P right, FLOAT_P right2, FLOAT_P dx_);

FLOAT_P upwind_first_derivative_y(FLOAT_P **array, FLOAT_P **velocity, int i, int j, FLOAT_P dy, int ny);
FLOAT_P upwind_first_derivative_z(FLOAT_P **array, FLOAT_P **velocity, int i, int j, FLOAT_P dz, int nz);

FLOAT_P central_first_derivative_y(FLOAT_P **array, int i, int j, FLOAT_P dy, int ny);
FLOAT_P central_first_derivative_z(FLOAT_P **array, int i, int j, FLOAT_P dz, int nz);

FLOAT_P central_second_derivative_y(FLOAT_P **array, int i, int j, FLOAT_P dy, int ny);
FLOAT_P central_second_derivative_z(FLOAT_P **array, int i, int j, FLOAT_P dz, int nz);

#endif // SPACIAL_DERIVATIVES_H__