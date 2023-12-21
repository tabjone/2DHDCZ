#include "global_float_precision.h"

void extrapolate_1D_array_constant_down(FLOAT_P *array, int nz_ghost)
{
    /*
    Extrapolates the ghost cells of a 1D array using constant extrapolation in the downward direction.

    Parameters
    ----------
    array : FLOAT_P*
        A pointer to the array to be extrapolated.
    grid_info : GridInfo
        A pointer to the GridInfo struct.
    */

    // Extrapolating
    for (int i = 0; i < nz_ghost; i++)
    {
        array[i] = array[nz_ghost];
    }
}