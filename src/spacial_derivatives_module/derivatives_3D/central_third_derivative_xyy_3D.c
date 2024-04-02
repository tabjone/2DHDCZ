#include "global_float_precision.h"
#include "periodic_boundary.h"

FLOAT_P central_third_derivative_xyy_3D(FLOAT_P ***array, int i, int j, int k, int nx, int ny, FLOAT_P one_over_8dxdydy)
{
    /*
    Calculates the central third derivative of a 3D array at a point (i, j, k) in the xyy direction.

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
    one_over_8dxdydy : FLOAT_P
        1/(8dxdydy).
    */

    int k_plus = periodic_boundary(k + 1, nx);
    int k_minus = periodic_boundary(k - 1, nx);
    int j_plus2 = periodic_boundary(j + 2, ny);
    int j_minus2 = periodic_boundary(j - 2, ny);

    return
    (
        array[i][j_plus2][k_plus]  - 2.0*array[i][j_plus2][k]  + array[i][j_plus2][k_minus]
      - array[i][j_minus2][k_plus] + 2.0*array[i][j_minus2][k] - array[i][j_minus2][k_minus]
    )* one_over_8dxdydy;
}