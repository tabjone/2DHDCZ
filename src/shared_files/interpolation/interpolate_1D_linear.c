#include "global_parameters.h"

FLOAT_P interpolate_1D_linear(FLOAT_P x0, FLOAT_P x1, FLOAT_P y0, FLOAT_P y1, FLOAT_P x) 
{
    /*
    Finds the value of y at x using linear interpolation on the interval [x0, x1] with corresponding values [y0, y1].

    Parameters
    ----------
    x0 : FLOAT_P
        The x value of the first point.
    x1 : FLOAT_P
        The x value of the second point.
    y0 : FLOAT_P
        The y value of the first point.
    y1 : FLOAT_P
        The y value of the second point.
    */
    return y0 + ((x - x0) * (y1 - y0) / (x1 - x0));
}