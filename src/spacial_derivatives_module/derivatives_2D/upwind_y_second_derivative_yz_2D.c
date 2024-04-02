#include "global_float_precision.h"
#include "global_parameters.h"
#include "periodic_boundary.h"

FLOAT_P upwind_y_second_derivative_yz_2D(FLOAT_P **array, FLOAT_P **velocity, int i, int j, int ny, FLOAT_P one_over_2dydz)
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
    ny : int
        The number of points in the y-direction.
    one_over_2dydz : FLOAT_P
        1/(2*dy*dz)
    */

    int j_minus = periodic_boundary(j - 1, ny);
    int j_plus = periodic_boundary(j + 1, ny);

    if (velocity[i][j] >= 0)
    {
        return (array[i+1][j] - array[i-1][j] - array[i+1][j_minus] + array[i-1][j_minus]) * one_over_2dydz;
    }
    else
    {
        return (array[i+1][j_plus] - array[i-1][j_plus] - array[i+1][j] + array[i-1][j]) * one_over_2dydz;
    }
}