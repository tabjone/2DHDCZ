#include "global_parameters.h"
#include "global_constants.h"
#include "global_float_precision.h"

#include <math.h>

FLOAT_P get_k_value(FLOAT_P r)
{
    return (tanh((r - 0.70*R_SUN)/(0.03*R_SUN)) + 1.0)/2.0 * SUPERAD_PARAM;
}
