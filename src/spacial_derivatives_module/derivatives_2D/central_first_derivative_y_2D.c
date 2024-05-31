#include "global_float_precision.h"
#include "periodic_boundary.h"
#include "global_parameters.h"

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

    #if CENTRAL_ORDER == 2
        int j_minus = periodic_boundary(j - 1, ny);
        int j_plus = periodic_boundary(j + 1, ny);
        return (array[i][j_plus] - array[i][j_minus]) * one_over_2dy;
    #elif CENTRAL_ORDER == 4
        int j_minus_2 = periodic_boundary(j - 2, ny);
        int j_minus_1 = periodic_boundary(j - 1, ny);
        int j_plus_1 = periodic_boundary(j + 1, ny);
        int j_plus_2 = periodic_boundary(j + 2, ny);
        return (array[i][j_minus_2] - 8.0 * array[i][j_minus_1] + 8.0 * array[i][j_plus_1] - array[i][j_plus_2]) * one_over_2dy / 6.0;
    #endif // CENTRAL_ORDER
}