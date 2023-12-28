#include "global_float_precision.h"
#include "global_parameters.h"

static inline int periodic_boundary(int i, int limit) {
    return (i + limit-1) % (limit-1);}

FLOAT_P upwind_first_derivative_x_3D(FLOAT_P ***array, FLOAT_P ***velocity, int i, int j, int k, int nx, FLOAT_P one_over_dx, FLOAT_P one_over_2dx)
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
            return (array[i][j][k] - array[i][j][k_minus]) * one_over_dx;
        }
        else
        {
            return (array[i][j][k_plus] - array[i][j][k]) * one_over_dx;
        }
    #elif UPWIND_ORDER == 2
        int k_minus2 = periodic_boundary(k-2, nx);
        int k_plus2 = periodic_boundary(k+2, nx);

        if (velocity[i][j][k] >= 0)
        {
            return (3.0*array[i][j][k] -4.0*array[i][j][k_minus] + array[i][j][k_minus2]) * one_over_2dx;
        }
        else
        {
            return (-3.0*array[i][j][k] + 4.0*array[i][j][k_plus] - array[i][j][k_plus2]) * one_over_2dx;
        }
    #endif // UPWIND_ORDER
}