double backward_first_derivative_second_order(double centre, double left, double left2, double dx_)
{
    /*
    Calculates the first derivate of a function using the second order backward difference method.

    Parameters
    ----------
    left : double
        The value of f(x-dx)
    left2: double
        The value of f(x-2dx)
    centre : double
        The value of f(x)

    Returns
    -------
    double
        The first derivate of f(x)
    */
    return (3*centre - 4*left + left2) / (2 * dx_);
}