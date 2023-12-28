#include "global_float_precision.h"

static inline int periodic_boundary(int i, int limit) {
    return (i + limit-1) % (limit-1);}

FLOAT_P central_second_derivative_y_2D(FLOAT_P **array, int i, int j, int ny, FLOAT_P one_over_dydy)
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
    ny : int
        The number of points in the y-direction.
    one_over_dydy : FLOAT_P
        1/dy^2.
    */

    int j_minus = periodic_boundary(j - 1, ny);
    int j_plus = periodic_boundary(j + 1, ny);

    return (array[i][j_plus] - 2.0*array[i][j] + array[i][j_minus]) * one_over_dydy;
}