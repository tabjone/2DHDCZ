#include "extrapolation.h"

void extrapolate_1D_array_constant_up(FLOAT_P *array, int nz_full, int nz_ghost)
{
    /*
    Extrapolates the ghost cells of a 1D array using constant extrapolation in the upward direction.

    Parameters
    ----------
    array : FLOAT_P*
        A pointer to the array to be extrapolated.
    grid_info : GridInfo
        A pointer to the GridInfo struct.
    */

    // Extrapolating
    for (int i = nz_full-nz_ghost; i < nz_full; i++)
    {
        array[i] = array[nz_full-nz_ghost-1];
    }
}