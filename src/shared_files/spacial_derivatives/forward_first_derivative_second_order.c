#include "global_parameters.h"

FLOAT_P forward_first_derivative_second_order(FLOAT_P centre, FLOAT_P right, FLOAT_P right2, FLOAT_P dx_)
{
    /*
    Calculates the first derivate of a function using the second order forward difference method.

    Parameters
    ----------
    centre : FLOAT_P
        The value of f(x)
    right : FLOAT_P
        The value of f(x+dx)
    right2 : FLOAT_P
        The value of f(x+2*dx)
    dx : FLOAT_P
        The step size
    
    Returns
    -------
    FLOAT_P
        The first derivate of f(x)
    */
    return (-3*centre + 4*right - right2) / (2 * dx_);
}