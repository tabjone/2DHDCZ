#include "extrapolation.h"

#if DIMENSIONS == 1
void extrapolate_1D_array_constant_up(FLOAT_P *array, struct GridInfo *grid_info)
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

    // Getting grid info
    int nz_ghost = grid_info->nz_ghost;
    int nz_full = grid_info->nz_full;

    // Extrapolating
    for (int i = nz_full-nz_ghost; i < nz_full; i++)
    {
        array[i] = array[nz_full-nz_ghost-1];
    }
}
#endif // DIMENSIONS