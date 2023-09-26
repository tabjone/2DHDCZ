#include "global_parameters.h"
#include "global_constants.h"

FLOAT_P get_k_value(FLOAT_P r)
{
    if (r>=0.7*R_SUN)
    {
        return K;
    }
    else
    {
        return 0.0;
    }
}