#include "global_float_precision.h"
#include "periodic_boundary.h"

FLOAT_P central_third_derivative_xxy_3D(FLOAT_P ***array, int i, int j, int k, int nx, int ny, FLOAT_P one_over_8dxdxdy)
{
    /*
    Calculates the central third derivative of a 3D array at a point (i, j, k) in the xxy direction.

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
    one_over_8dxdxdy : FLOAT_P
        1/(8dxdxdy).
    */

    int k_plus2 = periodic_boundary(k + 2, nx);
    int k_minus2 = periodic_boundary(k - 2, nx);
    int j_plus = periodic_boundary(j + 1, ny);
    int j_minus = periodic_boundary(j - 1, ny);

    return
    (
        array[i][j_plus][k_plus2]  - 2.0*array[i][j_plus][k]  + array[i][j_plus][k_minus2]
      - array[i][j_minus][k_plus2] + 2.0*array[i][j_minus][k] - array[i][j_minus][k_minus2]
    )* one_over_8dxdxdy;
}