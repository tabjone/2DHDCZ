#include "global_parameters.h"


FLOAT_P backward_first_derivative_second_order(FLOAT_P centre, FLOAT_P left, FLOAT_P left2, FLOAT_P dx_)
{
    /*
    Calculates the first derivate of a function using the second order backward difference method.

    Parameters
    ----------
    left : FLOAT_P
        The value of f(x-dx)
    left2: FLOAT_P
        The value of f(x-2dx)
    centre : FLOAT_P
        The value of f(x)

    Returns
    -------
    FLOAT_P
        The first derivate of f(x)
    */
    return (3*centre - 4*left + left2) / (2 * dx_);
}