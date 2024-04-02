#include "global_float_precision.h"
#include "periodic_boundary.h"

FLOAT_P central_third_derivative_y_2D(FLOAT_P **array, int i, int j, int ny, FLOAT_P one_over_2dydydy)
{
    /*
    Calculates the central third derivative of a 2D array at a point (i, j) in the y direction.

    Parameters
    ----------
    array : FLOAT_P**
        A pointer to the 2D array.
    i : int
        The i index of the point.
    j : int
        The j index of the point.
    ny : int
        The number of points in the y-direction.
    one_over_2dydydy : FLOAT_P
        1/(2dydydy).
    */

    int j_plus = periodic_boundary(j + 1, ny);
    int j_minus = periodic_boundary(j - 1, ny);
    int j_plus2 = periodic_boundary(j + 2, ny);
    int j_minus2 = periodic_boundary(j - 2, ny);

    return (array[i][j_plus2]-array[i][j_minus2] + 2.0*array[i][j_minus] - 2.0*array[i][j_plus]) * one_over_2dydydy;
}