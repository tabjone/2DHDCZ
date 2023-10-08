#include "extrapolation.h"

void extrapolate_3D_array_constant_down(FLOAT_P ***array, struct GridInfo *grid_info)
{
    /*
    Extrapolates the ghost cells of a 3D array using constant extrapolation in the downward direction.

    Parameters
    ----------
    array : FLOAT_P***
        A pointer to the array to be extrapolated.
    grid_info : GridInfo
        A pointer to the GridInfo struct.
    */

    // Getting grid info
    int nx = grid_info->nx;
    int ny = grid_info->ny;
    int nz_ghost = grid_info->nz_ghost;

    // Extrapolating
    for (int i = 0; i < nz_ghost; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                array[i][j][k] = array[nz_ghost][j][k];
            }
        }
    }
}