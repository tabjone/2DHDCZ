#include "spacial_derivatives_3D.h"

FLOAT_P central_second_derivative_x_3D(FLOAT_P ***array, int i, int j, int k, FLOAT_P dx, int nx)
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

    return (array[i][j][k_plus] - 2.0*array[i][j][k] + array[i][j][k_minus]) / (dx * dx);
}