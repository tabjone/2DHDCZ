#include "spacial_derivatives_3D.h"

FLOAT_P upwind_first_derivative_z_3D(FLOAT_P ***array, FLOAT_P ***velocity, int i, int j, int k, FLOAT_P dz, int nz)
{
    /*
    Calculates the central first derivative of a 3D array at a point (i, j, k) in the z-direction.

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
    dz : FLOAT_P
        The spacing between the points.
    nz : int
        The number of points in the z-direction.
    */

    #if UPWIND_ORDER == 1
        if (velocity[i][j][k] >= 0)
        {
            return (array[i][j][k] - array[i-1][j][k])/dz;
        }
        else
        {
            return (array[i+1][j][k] - array[i][j][k])/dz;
        }
    #elif UPWIND_ORDER == 2
        if (velocity[i][j][k] >= 0)
        {
            return (3.0*array[i][j][k] -4.0*array[i-1][j][k] + array[i-2][j][k])/(2.0*dz);
        }
        else
        {
            return (-3.0*array[i][j][k] + 4.0*array[i+1][j][k] - array[i+2][j][k])/(2.0*dz);
        }
    #endif // UPWIND_ORDER
}