#include "spacial_derivatives.h"

FLOAT_P central_second_derivative_yz(FLOAT_P **array, struct PrecalculatedVariables *precalc, int i, int j)
{
    /*
    Calculates the central second derivative of a 2D array at a point (i, j) in the yz-direction.

    Parameters
    ----------
    array : FLOAT_P**
        A pointer to the 2D array.
    i : int
        The i index of the point.
    j : int
        The j index of the point.
    dy : FLOAT_P
        The spacing between the points.
    ny : int
        The number of points in the y-direction.
    */

    int j_minus = precalc->j_minus[j];
    int j_plus = precalc->j_plus[j];

    return (array[i+1][j_plus] - array[i+1][j_minus] - array[i-1][j_plus] + array[i-1][j_minus]) * precalc->one_over_4dydz;
}