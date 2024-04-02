#include "global_float_precision.h"
#include "periodic_boundary.h"

FLOAT_P central_second_derivative_xz_3D(FLOAT_P ***array, int i, int j, int k, int nx, FLOAT_P one_over_4dxdz)
{
    /*
    Calculates the central second derivative of a 3D array at a point (i, j, k) in the xz-direction.

    Parameters
    ----------
    array : FLOAT_P***
        A pointer to the 3D array.
    i : int
        The i index of the point.
    j : int
        The j index of the point.
    nx : int
        The number of points in the x-direction.
    one_over_4dxdz : FLOAT_P
        1/(4dxdz).
    */

    int k_minus = periodic_boundary(k-1, nx);
    int k_plus = periodic_boundary(k+1, nx);

    return (array[i+1][j][k_plus] - array[i+1][j][k_minus] - array[i-1][j][k_plus] + array[i-1][j][k_minus]) * one_over_4dxdz;
}