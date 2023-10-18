#include "global_parameters.h"
#include "global_constants.h"

FLOAT_P get_k_value(FLOAT_P r)
{
    if (r>=0.7*R_SUN)
    {
        return K;
    }
    else if (0.65*R_SUN<r<0.7*R_SUN)
    {
        return K*(r-0.65*R_SUN)/(0.7*R_SUN-0.65*R_SUN);
    }
    else
    {
        return 0.0;
    }
}