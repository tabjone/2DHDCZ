#include "spacial_derivatives.h"

FLOAT_P central_third_derivative_yzz(FLOAT_P **array, struct PrecalculatedVariables *precalc, int i, int j)
{
    /*
    Calculates the central third derivative of a 2D array at a point (i, j) in the yzz-direction.

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

    int j_plus = precalc->j_plus[j];
    int j_minus = precalc->j_minus[j];

    return (array[i+2][j_plus] - 2.0*array[i][j_plus] +array[i-2][j_plus] - array[i+2][j_minus]  + 2.0*array[i][j_minus] - array[i-2][j_minus]) * precalc->one_over_8dydzdz;
}