#include "global_float_precision.h"

FLOAT_P central_third_derivative_z_2D(FLOAT_P **array, int i, int j, FLOAT_P one_over_2dzdzdz)
{
    /*
    Calculates the central third derivative of a 2D array at a point (i, j) in the z direction.

    Parameters
    ----------
    array : FLOAT_P**
        A pointer to the 2D array.
    i : int
        The i index of the point.
    j : int
        The j index of the point.
    one_over_2dzdzdz : FLOAT_P
        1/(2dzdzdz).
    */

    return (array[i+2][j]-array[i-2][j] + 2.0*array[i-1][j] - 2.0*array[i+1][j]) * one_over_2dzdzdz;
}