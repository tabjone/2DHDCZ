#include "global_parameters.h"
#include "global_constants.h"

double get_k_value(double r)
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