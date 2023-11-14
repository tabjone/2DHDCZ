#include "spacial_derivatives.h"

FLOAT_P central_third_derivative_yyz(FLOAT_P **array, struct PrecalculatedVariables *precalc, int i, int j)
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

    int j_plus2 = precalc->j_plus2[j];
    int j_minus2 = precalc->j_minus2[j];

    return (array[i+1][j_plus2]-2.0*array[i+1][j]+array[i+1][j_minus2] - array[i-1][j_plus2]+2.0*array[i-1][j]-array[i-1][j_minus2]) * precalc->one_over_8dydydz;
}