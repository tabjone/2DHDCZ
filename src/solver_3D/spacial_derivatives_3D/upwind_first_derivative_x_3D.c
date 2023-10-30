#include "spacial_derivatives_3D.h"

FLOAT_P upwind_first_derivative_x_3D(FLOAT_P ***array, FLOAT_P ***velocity, int i, int j, int k, FLOAT_P dx, int nx)
{
    /*
    Calculates the central first derivative of a 2D array at a point (i, j) in the x-direction.

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
    dx : FLOAT_P
        The spacing between the points.
    nx : int
        The number of points in the x-direction.
    */
    int k_minus = periodic_boundary(k-1, nx);
    int k_plus = periodic_boundary(k+1, nx);

    #if UPWIND_ORDER == 1
        if (velocity[i][j][k] >= 0)
        {
            return (array[i][j][k] - array[i][j][k_minus])/dx;
        }
        else
        {
            return (array[i][j][k_plus] - array[i][j][k])/dx;
        }
    #elif UPWIND_ORDER == 2
        int k_minus2 = periodic_boundary(k-2, nx);
        int k_plus2 = periodic_boundary(k+2, nx);

        if (velocity[i][j][k] >= 0)
        {
            return (3.0*array[i][j][k] -4.0*array[i][j][k_minus] + array[i][j][k_minus2])/(2.0*dx);
        }
        else
        {
            return (-3.0*array[i][j][k] + 4.0*array[i][j][k_plus] - array[i][j][k_plus2])/(2.0*dx);
        }
    #endif // UPWIND_ORDER
}