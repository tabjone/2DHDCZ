#include "global_float_precision.h"

FLOAT_P central_second_derivative_z_1D(FLOAT_P *array, int i, FLOAT_P one_over_dzdz)
{
    /*
    Calculates the central second derivative of a 1D array at a point (i) in the z-direction.

    Parameters
    ----------
    array : FLOAT_P*
        A pointer to the 1D array.
    i : int
        The i index of the point.
    one_over_dzdz : FLOAT_P
        1/dz^2.
    */
    return (array[i+1] - 2.0*array[i] + array[i-1]) * one_over_dzdz;
}