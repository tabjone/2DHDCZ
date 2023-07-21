double backward_first_derivative_first_order(double centre, double left, double dx)
{
    /*
    Calculates the first derivate of a function using the first order backward difference method.

    Parameters
    ----------
    left : double
        The value of f(x-dx)
    centre : double
        The value of f(x)

    Returns
    -------
    double
        The first derivate of f(x)
    */
    return (centre - left) / dx;
}