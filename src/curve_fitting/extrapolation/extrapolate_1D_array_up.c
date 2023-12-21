#include "extrapolation.h"
#include "global_float_precision.h"
#include "global_boundary.h"

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

    #if GHOST_CELLS_EXTRAPOLATION_VERTICAL == 0
        // Constant extrapolation
        extrapolate_1D_array_constant_up(array, nz_full, nz_ghost);
    #elif GHOST_CELLS_EXTRAPOLATION_VERTICAL == 1
        // Antisymmetric extrapolation
        extrapolate_1D_array_antisymmetric_up(array, nz_full, nz_ghost);
    #endif // EXTRAPOLATE_GHOST_CELLS
}