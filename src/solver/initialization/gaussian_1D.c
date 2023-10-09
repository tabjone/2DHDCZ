#include "initialization.h"

FLOAT_P gaussian_1D(FLOAT_P x, FLOAT_P x0, FLOAT_P sigma_x, FLOAT_P A) 
{
    /*
    Returns the value of a 1D Gaussian function at a given point.

    Parameters
    ----------
    x : FLOAT_P
        The x-coordinate of the point.
    x0 : FLOAT_P
        The x-coordinate of the center of the Gaussian.
    sigma_x : FLOAT_P
        The standard deviation of the Gaussian in the x-direction.
    A : FLOAT_P
        The amplitude of the Gaussian.

    Returns
    -------
    FLOAT_P
        The value of the Gaussian at the given point.
    */

    return A * exp(-(x - x0) * (x - x0) / (2 * sigma_x * sigma_x));
}