#include "global_float_precision.h"

static inline int periodic_boundary(int i, int limit) {
    return (i + limit-1) % (limit-1);}

void central_third_derivative_xxz_3D(FLOAT_P ***array, int i, int j, int k, int nx, FLOAT_P one_over_8dxdxdz)
{
    /*
    Calculates the central third derivative of a 3D array at a point (i, j, k) in the xxz direction.

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
    one_over_8dxdxdz : FLOAT_P
        1/(8dxdxdz).
    */

    int k_plus2 = periodic_boundary(k + 2, nx);
    int k_minus2 = periodic_boundary(k - 2, nx);

    return
    (
        array[i+1][j][k_plus2] - 2.0*array[i+1][j][k] + array[i+1][j][k_minus2]
      - array[i-1][j][k_plus2] + 2.0*array[i-1][j][k] - array[i-1][j][k_minus2]
    )* one_over_8dxdxdz;
}