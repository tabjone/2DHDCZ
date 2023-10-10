#include "extrapolation.h"

#if DIMENSIONS == 1
void extrapolate_1D_array_down(FLOAT_P *array, struct GridInfo *grid_info)
{
    /*
    Extrapolates the ghost cells of a 1D array in the downward direction.

    Parameters
    ----------
    array : FLOAT_P*
        A pointer to the array to be extrapolated.
    grid_info : GridInfo
        A pointer to the GridInfo struct.
    */

    #if EXTRAPOLATE_GHOST_CELLS == 0
        // Constant extrapolation
        extrapolate_1D_array_constant_down(array, grid_info);
    #else
        #error "Only constant extrapolation implemented."
    #endif // EXTRAPOLATE_GHOST_CELLS == 0
}
#endif // DIMENSIONS == 1