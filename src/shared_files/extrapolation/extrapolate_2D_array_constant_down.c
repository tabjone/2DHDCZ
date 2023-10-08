#include "extrapolation.h"

void extrapolate_2D_array_constant_down(FLOAT_P **array, struct GridInfo *grid_info)
{
    /*
    Extrapolates the ghost cells of a 2D array using constant extrapolation in the downward direction.

    Parameters
    ----------
    array : FLOAT_P**
        A pointer to the array to be extrapolated.
    grid_info : GridInfo
        A pointer to the GridInfo struct.
    */

    // Getting grid info
    int ny = grid_info->ny;
    int nz_ghost = grid_info->nz_ghost;

    // Extrapolating
    for (int i = 0; i < nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            array[i][j] = array[nz_ghost][j];
        }
    }
}