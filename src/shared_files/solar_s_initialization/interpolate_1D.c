double interpolate_1D(double x0, double x1, double y0, double y1, double x) 
{
    /*
    Finds the value of y at x using linear interpolation on the interval [x0, x1] with corresponding values [y0, y1].

    Parameters
    ----------
    x0 : double
        The x value of the first point.
    x1 : double
        The x value of the second point.
    y0 : double
        The y value of the first point.
    y1 : double
        The y value of the second point.
    */
    return y0 + ((x - x0) * (y1 - y0) / (x1 - x0));
}