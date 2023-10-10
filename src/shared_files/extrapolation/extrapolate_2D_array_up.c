#include "extrapolation.h"

#if DIMENSIONS == 2
void extrapolate_2D_array_up(FLOAT_P **array, struct GridInfo *grid_info)
{
    /*
    Extrapolates the ghost cells of a 2D array in the upward direction.

    Parameters
    ----------
    array : FLOAT_P**
        A pointer to the array to be extrapolated.
    grid_info : GridInfo
        A pointer to the GridInfo struct.
    */

    #if EXTRAPOLATE_GHOST_CELLS == 0
        // Constant extrapolation
        extrapolate_2D_array_constant_up(array, grid_info);
    #else
        #error "Only constant extrapolation implemented."
    #endif // EXTRAPOLATE_GHOST_CELLS == 0
}
#endif // DIMENSIONS == 2