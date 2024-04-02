#include "global_float_precision.h"
#include "periodic_boundary.h"

FLOAT_P central_third_derivative_xyz_3D(FLOAT_P ***array, int i, int j, int k, int nx, int ny, FLOAT_P one_over_8dxdydz)
{
    /*
    Calculates the central third derivative of a 3D array at a point (i, j, k) in the xyz direction.

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
    one_over_8dxdydz : FLOAT_P
        1/(8dxdydz).
    */

    int j_plus = periodic_boundary(j + 1, ny);
    int j_minus = periodic_boundary(j - 1, ny);
    int k_plus = periodic_boundary(k + 1, nx);
    int k_minus = periodic_boundary(k - 1, nx);

    return
    (
        array[i+1][j_plus][k_plus]  - array[i+1][j_plus][k_minus]
      - array[i+1][j_minus][k_plus] + array[i+1][j_minus][k_minus]
      - array[i-1][j_plus][k_plus]  + array[i-1][j_plus][k_minus]
      + array[i-1][j_minus][k_plus] - array[i-1][j_minus][k_minus]
    )* one_over_8dxdydz;
}