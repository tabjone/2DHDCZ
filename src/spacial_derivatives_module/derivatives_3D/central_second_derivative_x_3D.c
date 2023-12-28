#include "global_float_precision.h"

static inline int periodic_boundary(int i, int limit) {
    return (i + limit-1) % (limit-1);}

FLOAT_P central_second_derivative_x_3D(FLOAT_P ***array, int i, int j, int k, int nx, FLOAT_P one_over_dxdx)
{
    /*
    Calculates the central second derivative of a 2D array at a point (i, j) in the y-direction.

    Parameters
    ----------
    array : FLOAT_P**
        A pointer to the 2D array.
    i : int
        The i index of the point.
    j : int
        The j index of the point.
    dy : FLOAT_P
        The spacing between the points.
    ny : int
        The number of points in the y-direction.
    */

    int k_minus = periodic_boundary(k-1, nx);
    int k_plus = periodic_boundary(k+1, nx);

    return (array[i][j][k_plus] - 2.0*array[i][j][k] + array[i][j][k_minus]) * one_over_dxdx;
}