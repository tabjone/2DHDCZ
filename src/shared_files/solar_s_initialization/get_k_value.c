#include "global_parameters.h"
#include "global_constants.h"

FLOAT_P get_k_value(FLOAT_P r)
{
    FLOAT_P return_value = 0.0;
    if (r>=CZ_START*R_SUN)
    {
        return_value = SUPERAD_PARAM;
    }
    else if (r < CZ_START*R_SUN && r > 0.65*R_SUN)
    {
        return return_value = SUPERAD_PARAM*(r-0.65*R_SUN)/(0.7*R_SUN-0.65*R_SUN);
    }
    else
    {
        return_value = 0.0;
    }
    return return_value;
}
