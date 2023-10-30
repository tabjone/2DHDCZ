#include "spacial_derivatives_3D.h"

FLOAT_P upwind_first_derivative_y_3D(FLOAT_P ***array, FLOAT_P ***velocity, int i, int j, int k, FLOAT_P dy, int ny)
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
    int j_plus = periodic_boundary(j+1, ny);

    #if UPWIND_ORDER == 1
        if (velocity[i][j][k] >= 0)
        {
            return (array[i][j][k] - array[i][j_minus][k])/dy;
        }
        else
        {
            return (array[i][j_plus][k] - array[i][j][k])/dy;
        }
    #elif UPWIND_ORDER == 2
        int j_minus2 = periodic_boundary(j-2, ny);
        int j_plus2 = periodic_boundary(j+2, ny);

        if (velocity[i][j][k] >= 0)
        {
            return (3.0*array[i][j][k] -4.0*array[i][j_minus][k] + array[i][j_minus2][k])/(2.0*dy);
        }
        else
        {
            return (-3.0*array[i][j][k] + 4.0*array[i][j_plus][k] - array[i][j_plus2][k])/(2.0*dy);
        }
    #endif // UPWIND_ORDER
}