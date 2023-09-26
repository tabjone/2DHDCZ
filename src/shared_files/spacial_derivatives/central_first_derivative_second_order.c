#include "global_parameters.h"


FLOAT_P central_first_derivative_second_order(FLOAT_P left, FLOAT_P right, FLOAT_P dx_)
{
    /*
    Calculates the first derivate of a function using the second order central difference method.

    Parameters
    ----------
    left : FLOAT_P
        The value of f(x-dx)
    right : FLOAT_P
        The value of f(x+dx)
    dx : FLOAT_P
        The step size
    
    Returns
    -------
    FLOAT_P
        The first derivate of f(x)
    */
    return (right - left) / (2 * dx_);
}