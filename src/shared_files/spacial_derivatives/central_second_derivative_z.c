#include "spacial_derivatives.h"

FLOAT_P central_second_derivative_z(FLOAT_P **array, int i, int j, FLOAT_P dz, int nz)
{
    /*
    Calculates the central second derivative of a 2D array at a point (i, j) in the z-direction.

    Parameters
    ----------
    array : FLOAT_P**
        A pointer to the 2D array.
    i : int
        The i index of the point.
    j : int
        The j index of the point.
    dz : FLOAT_P
        The spacing between the points.
    nz : int
        The number of points in the z-direction.
    */

   #if CENTRAL_ORDER == 2
        return central_second_derivative_second_order(array[i][j], array[i-1][j], array[i+1][j], dz);
    #endif // CENTRAL_ORDER
}