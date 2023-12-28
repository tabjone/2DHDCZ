#include "global_float_precision.h"
#include "global_parameters.h"

FLOAT_P upwind_first_derivative_z_1D(FLOAT_P *array, FLOAT_P *velocity, int i, FLOAT_P one_over_dz, FLOAT_P one_over_2dz)
{
    /*
    Calculates the central first derivative of a 1D array at a point (i) in the z-direction.

    Parameters
    ----------
    array : FLOAT_P*
        A pointer to the 1D array.
    velocity: FLOAT_P*
        A pointer to the 1D array of velocity multiplied with this derivative.
    i : int
        The i index of the point.
    one_over_dz : FLOAT_P
        1/dz.
    one_over_2dz : FLOAT_P
        1/2dz.
    */

    #if UPWIND_ORDER == 1
        if (velocity[i] >= 0)
        {
            return (array[i] - array[i-1]) * one_over_dz;
        }
        else
        {
            return (array[i+1] - array[i]) * one_over_dz;
        }
    #elif UPWIND_ORDER == 2
        if (velocity[i] >= 0)
        {
            return (3.0*array[i] -4.0*array[i-1] + array[i-2]) * one_over_2dz;
        }
        else
        {
            return (-3.0*array[i] + 4.0*array[i+1] - array[i+2]) * one_over_2dz;
        }
    #endif // UPWIND_ORDER
}