#include "global_float_precision.h"
#include "periodic_boundary.h"

FLOAT_P central_third_derivative_yzz_3D(FLOAT_P ***array, int i, int j, int k, int ny, FLOAT_P one_over_8dydzdz)
{
    /*
    Calculates the central third derivative of a 3D array at a point (i, j, k) in the yzz direction.

    Parameters
    ----------
    array : FLOAT_P***
        A pointer to the 3D array.
    i : int
        The i index of the point.
    j : int
        The j index of the point.
    k : int
        The k index of the point.
    ny : int
        The number of points in the y-direction.
    one_over_8dydzdz : FLOAT_P
        1/(8dydzdz).
    */

    int j_plus = periodic_boundary(j + 1, ny);
    int j_minus = periodic_boundary(j - 1, ny);

    return
    (
        array[i+2][j_plus][k]  - 2.0*array[i+2][j][k] + array[i+2][j_minus][k]
      - array[i-2][j_plus][k]  + 2.0*array[i-2][j][k] - array[i-2][j_minus][k]
    )* one_over_8dydzdz;
}