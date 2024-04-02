#include "global_float_precision.h"
#include "periodic_boundary.h"

FLOAT_P central_third_derivative_xzz_3D(FLOAT_P ***array, int i, int j, int k, int nx, FLOAT_P one_over_8dxdzdz)
{
    /*
    Calculates the central third derivative of a 3D array at a point (i, j, k) in the xzz direction.

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
    nx : int
        The number of points in the x-direction.
    ny : int
        The number of points in the y-direction.
    one_over_8dxdzdz : FLOAT_P
        1/(8dxdzdz).
    */

    int k_plus = periodic_boundary(k + 1, nx);
    int k_minus = periodic_boundary(k - 1, nx);

    return
    (
        array[i+2][j][k_plus] - 2.0*array[i+2][j][k] + array[i+2][j][k_minus]
      - array[i-2][j][k_plus] + 2.0*array[i-2][j][k] - array[i-2][j][k_minus]
    )* one_over_8dxdzdz;
}