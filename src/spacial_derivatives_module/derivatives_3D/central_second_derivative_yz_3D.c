#include "global_float_precision.h"

static inline int periodic_boundary(int i, int limit) {
    return (i + limit-1) % (limit-1);}

FLOAT_P central_second_derivative_yz_3D(FLOAT_P ***array, int i, int j, int k, int ny, FLOAT_P one_over_4dydz)
{
    /*
    Calculates the central second derivative of a 3D array at a point (i, j, k) in the yz-direction.

    Parameters
    ----------
    array : FLOAT_P***
        A pointer to the 3D array.
    i : int
        The i index of the point.
    j : int
        The j index of the point.
    ny : int
        The number of points in the y-direction.
    one_over_4dydz : FLOAT_P
        1/(4dydz).
    */

    int j_minus = periodic_boundary(j-1, ny);
    int j_plus = periodic_boundary(j+1, ny);

    return (array[i+1][j_plus][k] - array[i+1][j_minus][k] - array[i-1][j_plus][k] + array[i-1][j_minus][k]) * one_over_4dydz;
}