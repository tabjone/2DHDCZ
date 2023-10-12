#include "spacial_derivatives.h"

FLOAT_P upwind_first_derivative_y(FLOAT_P **array, FLOAT_P **velocity, int i, int j, FLOAT_P dy, int ny)
{
    /*
    Calculates the central first derivative of a 2D array at a point (i, j) in the y-direction.

    Parameters
    ----------
    array : FLOAT_P**
        A pointer to the 2D array.
    velocity: FLOAT_P**
        A pointer to the 2D array of velocity multiplied with this derivative.
    i : int
        The i index of the point.
    j : int
        The j index of the point.
    dy : FLOAT_P
        The spacing between the points.
    ny : int
        The number of points in the z-direction.
    */

    int j_minus = periodic_boundary(j-1, ny);
    int j_minus2 = periodic_boundary(j-2, ny);
    int j_plus = periodic_boundary(j+1, ny);
    int j_plus2 = periodic_boundary(j+2, ny);

    #if UPWIND_ORDER == 1
        if (velocity[i][j] >= 0)
        {
            return backward_first_derivative_first_order(array[i][j], array[i][j_minus], dy);
        }
        else
        {
            return forward_first_derivative_first_order(array[i][j], array[i][j_plus], dy);
        }
    #elif UPWIND_ORDER == 2
        if (velocity[i][j] >= 0)
        {
            return backward_first_derivative_second_order(array[i][j], array[i][j_minus], array[i][j_minus2], dy);
        }
        else
        {
            return forward_first_derivative_second_order(array[i][j], array[i][j_plus], array[i][j_plus2], dy);
        }
    #endif // UPWIND_ORDER
}