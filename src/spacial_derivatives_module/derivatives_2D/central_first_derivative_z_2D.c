#include "global_float_precision.h"
#include "global_parameters.h"

FLOAT_P central_first_derivative_z_2D(FLOAT_P **array, int i, int j, FLOAT_P one_over_2dz)
{
    /*
    Calculates the central first derivative of a 2D array at a point (i, j) in the z-direction.

    Parameters
    ----------
    array : FLOAT_P**
        A pointer to the 2D array.
    i : int
        The i index of the point.
    j : int
        The j index of the point.
    one_over_2dz : FLOAT_P
        1/2dz.
    */
    #if CENTRAL_ORDER == 2
        return (array[i+1][j] - array[i-1][j]) * one_over_2dz;
    #elif CENTRAL_ORDER == 4
        return (array[i-2][j] - 8.0 * array[i-1][j] + 8.0 * array[i+1][j] - array[i+2][j]) * one_over_2dz / 6.0;
    #endif // CENTRAL_ORDER
}