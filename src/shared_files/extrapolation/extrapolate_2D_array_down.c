#include "extrapolation.h"

void extrapolate_2D_array_down(FLOAT_P **array, int nz_ghost, int ny)
{
    /*
    Extrapolates the ghost cells of a 2D array in the downward direction.

    Parameters
    ----------
    array : FLOAT_P**
        A pointer to the array to be extrapolated.
    grid_info : GridInfo
        A pointer to the GridInfo struct.
    */
    #if GHOST_CELLS_EXTRAPOLATION_VERTICAL == 0
        // Constant extrapolation
        extrapolate_2D_array_constant_down(array, nz_ghost, ny);
    #elif GHOST_CELLS_EXTRAPOLATION_VERTICAL == 1
        // Antisymmetric extrapolation
        extrapolate_2D_array_antisymmetric_down(array, nz_ghost, ny);
    #else
        #error "Extrapolating bla bla."
    #endif // EXTRAPOLATE_GHOST_CELLS
}