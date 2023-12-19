#include "derivatives_2D.h"

static inline int periodic_boundary(int i, int limit) {
    return (i + limit-1) % (limit-1);}

FLOAT_P central_first_derivative_y_2D(FLOAT_P **array, int i, int j, int ny, FLOAT_P one_over_2dy)
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
    ny : int
        The number of points in the y-direction.
    one_over_2dy : FLOAT_P
        1/2dy.

    */
    int j_minus = periodic_boundary(j - 1, ny);
    int j_plus = periodic_boundary(j + 1, ny);

    return (array[i][j_plus] - array[i][j_minus]) * one_over_2dy;
}