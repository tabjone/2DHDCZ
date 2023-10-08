#include "extrapolation.h"

void extrapolate_3D_array_constant_up(FLOAT_P ***array, struct GridInfo *grid_info)
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

    // Getting grid info
    int nx = grid_info->nx;
    int ny = grid_info->ny;
    int nz_ghost = grid_info->nz_ghost;
    int nz_full = grid_info->nz_full;

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