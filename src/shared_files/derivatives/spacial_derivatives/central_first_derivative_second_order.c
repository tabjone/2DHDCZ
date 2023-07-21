double central_first_derivative_second_order(double left, double right, double dx)
{
    /*
    Calculates the first derivate of a function using the second order central difference method.

    Parameters
    ----------
    left : double
        The value of f(x-dx)
    right : double
        The value of f(x+dx)
    dx : double
        The step size
    
    Returns
    -------
    double
        The first derivate of f(x)
    */
    return (right - left) / (2 * dx);
}