#include "global_parameters.h"


FLOAT_P central_second_derivative_second_order(FLOAT_P centre, FLOAT_P left, FLOAT_P right, FLOAT_P dx_)
{
    /*
    Calculates the second derivate of a function using the second order central difference method.

    Parameters
    ----------
    left : FLOAT_P
        The value of f(x-dx)
    right : FLOAT_P
        The value of f(x+dx)
    centre : FLOAT_P
        The value of f(x)
    dx : FLOAT_P
        The step size
    
    Returns
    -------
    FLOAT_P
        The second derivate of f(x)
    */
    return (right - 2*centre + left) / (dx_ * dx_);
}