double forward_first_derivative_first_order(double centre, double right, double dx)
{
    /*
    Calculates the first derivate of a function using the first order forward difference method.

    Parameters
    ----------
    right : double
        The value of f(x+dx)
    centre : double
        The value of f(x)

    Returns
    -------
    double
        The first derivate of f(x)
    */
    return (right - centre) / dx;
}