#ifndef SPACIAL_DERIVATIVES_H__
#define SPACIAL_DERIVATIVES_H__

#include "global_parameters.h"
#include "shared_files.h"
#include "../functions.h"

FLOAT_P upwind_first_derivative_y(FLOAT_P **array, FLOAT_P **velocity, struct PrecalculatedVariables *precalc, int i, int j);
FLOAT_P upwind_first_derivative_z(FLOAT_P **array, FLOAT_P **velocity, struct PrecalculatedVariables *precalc, int i, int j);

FLOAT_P central_first_derivative_y(FLOAT_P **array, struct PrecalculatedVariables *precalc, int i, int j);
FLOAT_P central_first_derivative_z(FLOAT_P **array, struct PrecalculatedVariables *precalc, int i, int j);

FLOAT_P central_second_derivative_y(FLOAT_P **array, struct PrecalculatedVariables *precalc, int i, int j);
FLOAT_P central_second_derivative_z(FLOAT_P **array, struct PrecalculatedVariables *precalc, int i, int j);

FLOAT_P central_second_derivative_yz(FLOAT_P **array, struct PrecalculatedVariables *precalc, int i, int j);

FLOAT_P central_third_derivative_yyz(FLOAT_P **array, struct PrecalculatedVariables *precalc, int i, int j);
FLOAT_P central_third_derivative_yzz(FLOAT_P **array, struct PrecalculatedVariables *precalc, int i, int j);



#endif // SPACIAL_DERIVATIVES_H__