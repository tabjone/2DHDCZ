#include "spacial_derivatives.h"

FLOAT_P upwind_first_derivative_y(FLOAT_P **array, FLOAT_P **velocity, struct PrecalculatedVariables *precalc, int i, int j)
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

    int j_minus = precalc->j_minus[j];
    int j_plus = precalc->j_plus[j];

    #if UPWIND_ORDER == 1
        if (velocity[i][j] >= 0)
        {
            return (array[i][j] - array[i][j_minus]) * precalc->one_over_dy;
        }
        else
        {
            return (array[i][j_plus] - array[i][j]) * precalc->one_over_dy;
        }
    #elif UPWIND_ORDER == 2
        int j_minus2 = precalc->j_minus2[j];
        int j_plus2 = precalc->j_plus2[j];

        if (velocity[i][j] >= 0)
        {
            return (3.0*array[i][j] -4.0*array[i][j_minus] + array[i][j_minus2]) * precalc->one_over_2dy;
        }
        else
        {
            return (-3.0*array[i][j] + 4.0*array[i][j_plus] - array[i][j_plus2]) * precalc->one_over_2dy;
        }
    #endif // UPWIND_ORDER
}