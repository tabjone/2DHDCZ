#include "global_float_precision.h"
#include "periodic_boundary.h"

FLOAT_P central_third_derivative_yzz_2D(FLOAT_P **array, int i, int j, int ny, FLOAT_P one_over_8dydzdz)
{
    /*
    Calculates the central third derivative of a 2D array at a point (i, j) in the yzz-direction.

    Parameters
    ----------
    array : FLOAT_P**
        A pointer to the 2D array.
    i : int
        The i index of the point.
    j : int
        The j index of the point.
    ny : int
        The number of points in the y-direction.
    one_over_8dydzdz : FLOAT_P
        1/(8dydzdz).
    */

    int j_plus = periodic_boundary(j + 1, ny);
    int j_minus = periodic_boundary(j - 1, ny);

    return (array[i+2][j_plus] - 2.0*array[i][j_plus] +array[i-2][j_plus] - array[i+2][j_minus]  + 2.0*array[i][j_minus] - array[i-2][j_minus]) * one_over_8dydzdz;
}