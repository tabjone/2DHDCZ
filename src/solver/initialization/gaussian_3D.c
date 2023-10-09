#include "initialization.h"

FLOAT_P gaussian_3D(FLOAT_P x, FLOAT_P y, FLOAT_P z, FLOAT_P x0, FLOAT_P y0, FLOAT_P z0, FLOAT_P sigma_x, FLOAT_P sigma_y, FLOAT_P sigma_z, FLOAT_P A) 
{
    /*
    Returns the value of a 3D Gaussian function at a given point.

    Parameters
    ----------
    x : FLOAT_P
        The x-coordinate of the point.
    y : FLOAT_P
        The y-coordinate of the point.
    z : FLOAT_P
        The z-coordinate of the point.
    x0 : FLOAT_P
        The x-coordinate of the center of the Gaussian.
    y0 : FLOAT_P
        The y-coordinate of the center of the Gaussian.
    z0 : FLOAT_P
        The z-coordinate of the center of the Gaussian.
    sigma_x : FLOAT_P
        The standard deviation of the Gaussian in the x-direction.
    sigma_y : FLOAT_P
        The standard deviation of the Gaussian in the y-direction.
    sigma_z : FLOAT_P
        The standard deviation of the Gaussian in the z-direction.
    A : FLOAT_P
        The amplitude of the Gaussian.

    Returns
    -------
    FLOAT_P
        The value of the Gaussian at the given point.
    */

    return A * exp(-(x - x0) * (x - x0) / (2 * sigma_x * sigma_x) 
                   -(y - y0) * (y - y0) / (2 * sigma_y * sigma_y)
                   -(z - z0) * (z - z0) / (2 * sigma_z * sigma_z));
}