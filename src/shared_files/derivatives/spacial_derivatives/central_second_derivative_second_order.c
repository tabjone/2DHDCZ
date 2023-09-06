double central_second_derivative_second_order(double centre, double left, double right, double dx_)
{
    /*
    Calculates the second derivate of a function using the second order central difference method.

    Parameters
    ----------
    left : double
        The value of f(x-dx)
    right : double
        The value of f(x+dx)
    centre : double
        The value of f(x)
    dx : double
        The step size
    
    Returns
    -------
    double
        The second derivate of f(x)
    */
    return (right - 2*centre + left) / (dx_ * dx_);
}