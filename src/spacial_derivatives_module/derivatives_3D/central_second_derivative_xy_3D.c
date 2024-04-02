#include "global_float_precision.h"
#include "periodic_boundary.h"

FLOAT_P central_second_derivative_xy_3D(FLOAT_P ***array, int i, int j, int k, int nx, int ny, FLOAT_P one_over_4dxdy)
{
    /*
    Calculates the central second derivative of a 3D array at a point (i, j, k) in the xy-direction.

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
    ny : int
        The number of points in the y-direction.
    one_over_4dxdy : FLOAT_P
        1/(4dxdy).
    */

    int k_minus = periodic_boundary(k-1, nx);
    int k_plus = periodic_boundary(k+1, nx);
    int j_minus = periodic_boundary(j-1, ny);
    int j_plus = periodic_boundary(j+1, ny);

    return (array[i][j_plus][k_plus] - array[i][j_plus][k_minus] - array[i][j_minus][k_plus] + array[i][j_minus][k_minus]) * one_over_4dxdy;
}