#include "global_float_precision.h"

void extrapolate_3D_array_constant_up(FLOAT_P ***array, int nz_full, int nz_ghost, int ny, int nx)
{
    /*
    Extrapolates the ghost cells of a 3D array using constant extrapolation in the upward direction.

    Parameters
    ----------
    array : FLOAT_P***
        A pointer to the array to be extrapolated.
    grid_info : GridInfo
        A pointer to the GridInfo struct.
    */

    // Extrapolating
    for (int i = nz_full-nz_ghost; i < nz_full; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                array[i][j][k] = array[nz_full-nz_ghost-1][j][k];
            }
        }
    }
}