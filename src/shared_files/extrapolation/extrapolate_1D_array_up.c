#include "extrapolation.h"

void extrapolate_1D_array_up(FLOAT_P *array, int nz_full, int nz_ghost)
{
    /*
    Extrapolates the ghost cells of a 1D array in the upward direction.

    Parameters
    ----------
    array : FLOAT_P*
        A pointer to the array to be extrapolated.
    grid_info : GridInfo
        A pointer to the GridInfo struct.
    */

    #if EXTRAPOLATE_GHOST_CELLS == 0
        // Constant extrapolation
        extrapolate_1D_array_constant_up(array, nz_full, nz_ghost);
    #else
        #error "Only constant extrapolation implemented."
    #endif // EXTRAPOLATE_GHOST_CELLS == 0
}