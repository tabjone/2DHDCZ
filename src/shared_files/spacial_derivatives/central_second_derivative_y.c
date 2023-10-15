#include "spacial_derivatives.h"

FLOAT_P central_second_derivative_y(FLOAT_P **array, int i, int j, FLOAT_P dy, int ny)
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

    int j_minus = periodic_boundary(j-1, ny);
    int j_plus = periodic_boundary(j+1, ny);

   #if CENTRAL_ORDER == 2
        return central_second_derivative_second_order(array[i][j], array[i][j_minus], array[i][j_plus], dy);
    #endif // CENTRAL_ORDER
}