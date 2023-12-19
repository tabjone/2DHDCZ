#include "derivatives_2D.h"

static inline int periodic_boundary(int i, int limit) {
    return (i + limit-1) % (limit-1);}

FLOAT_P central_third_derivative_yyz_2D(FLOAT_P **array, int i, int j, int ny, FLOAT_P one_over_8dydydz)
{
    /*
    Calculates the central third derivative of a 2D array at a point (i, j) in the yyz direction.

    Parameters
    ----------
    array : FLOAT_P**
        A pointer to the 2D array.
    i : int
        The i index of the point.
    j : int
        The j index of the point.
    ny : int
        The number of points in the y-direction.
    one_over_8dydydz : FLOAT_P
        1/(8dydydz).
    */

    int j_plus2 = periodic_boundary(j + 2, ny);
    int j_minus2 = periodic_boundary(j - 2, ny);

    return (array[i+1][j_plus2]-2.0*array[i+1][j]+array[i+1][j_minus2] - array[i-1][j_plus2]+2.0*array[i-1][j]-array[i-1][j_minus2]) * one_over_8dydydz;
}