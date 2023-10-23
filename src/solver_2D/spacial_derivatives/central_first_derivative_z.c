#include "spacial_derivatives.h"

FLOAT_P central_first_derivative_z(FLOAT_P **array, int i, int j, FLOAT_P dz, int nz)
{
    /*
    Calculates the central first derivative of a 2D array at a point (i, j) in the z-direction.

    Parameters
    ----------
    array : FLOAT_P**
        A pointer to the 2D array.
    i : int
        The i index of the point.
    j : int
        The j index of the point.
    dz : FLOAT_P
        The spacing between the points.
    nz : int
        The number of points in the z-direction.
    */
    return (array[i+1][j] - array[i-1][j]) / (2.0 * dz);
}