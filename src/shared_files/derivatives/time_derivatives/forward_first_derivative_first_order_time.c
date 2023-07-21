double forward_first_derivative_first_order_time(double current, double rhs, double dt)
{
    /*
    Calculates the value for the next timestep of a function using the first order Euler forward method.

    Parameters
    ----------
    current : double
        y_n
    rhs : double
        The value of f(x_n, t_n)
    dt : double
        The step size

    Returns
    -------
    double
        y_n+1
    */
    return current + dt * rhs;
}