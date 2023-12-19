#include "derivatives_2D.h"

FLOAT_P central_first_derivative_z_2D(FLOAT_P **array, int i, int j, FLOAT_P one_over_2dz)
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
    one_over_2dz : FLOAT_P
        1/2dz.
    */
    return (array[i+1][j] - array[i-1][j]) * one_over_2dz;
}