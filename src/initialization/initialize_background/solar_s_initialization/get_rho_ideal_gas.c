#include "global_float_precision.h"
#include "global_constants.h"

FLOAT_P get_rho_ideal_gas(FLOAT_P p0, FLOAT_P T0)
{
    return M_U * MU / K_B * p0/T0;
}