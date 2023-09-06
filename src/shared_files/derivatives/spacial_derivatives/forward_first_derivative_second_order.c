double forward_first_derivative_second_order(double centre, double right, double right2, double dx_)
{
    /*
    Calculates the first derivate of a function using the second order forward difference method.

    Parameters
    ----------
    centre : double
        The value of f(x)
    right : double
        The value of f(x+dx)
    right2 : double
        The value of f(x+2*dx)
    dx : double
        The step size
    
    Returns
    -------
    double
        The first derivate of f(x)
    */
    return (-3*centre + 4*right - right2) / (2 * dx_);
}