#include "global_parameters.h"


FLOAT_P backward_first_derivative_first_order(FLOAT_P centre, FLOAT_P left, FLOAT_P dx_)
{
    /*
    Calculates the first derivate of a function using the first order backward difference method.

    Parameters
    ----------
    left : FLOAT_P
        The value of f(x-dx)
    centre : FLOAT_P
        The value of f(x)

    Returns
    -------
    FLOAT_P
        The first derivate of f(x)
    */
    return (centre - left) / dx_;
}