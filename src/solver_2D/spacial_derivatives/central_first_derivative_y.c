#include "spacial_derivatives.h"

FLOAT_P central_first_derivative_y(FLOAT_P **array, int i, int j, FLOAT_P dy, int ny)
{
    /*
    Calculates the central first derivative of a 2D array at a point (i, j) in the y-direction.

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
    int j_minus = periodic_boundary(j-1, ny);
    int j_plus = periodic_boundary(j+1, ny);

    return (array[i][j_plus] - array[i][j_minus]) / (2.0 * dy);
}