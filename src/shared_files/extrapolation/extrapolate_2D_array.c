#include "global_parameters.h"

void extrapolate_2D_array(FLOAT_P **array, int nz_full, int nz_ghost, int ny)
{
    // Extrapolating p1 (constant extrapolation, move this to function later with possiblity for different extrapolation scheme)

    #if EXTAPOLATE_GHOST_CELLS == 0
    // Constant extrapolation
    for (int k = 0; k < nz_ghost; k++)
    {
        for (int j = 0; j < ny; j++)
        {
            array[k][j] = array[nz_ghost][j];
            array[nz_full-k-1][j] = array[nz_full-nz_ghost-1][j];
        }
    }
    #elif EXTAPOLATE_GHOST_CELLS == 1
    // Linear extrapolation (chatGPT)
    for (int k = 0; k < nz_ghost; k++)
    {
        for (int j = 0; j < ny; j++)
        {
            array[k][j] = 2*array[nz_ghost][j] - array[nz_ghost+1][j];
            array[nz_full-k-1][j] = 2*array[nz_full-nz_ghost-1][j] - array[nz_full-nz_ghost-2][j];
        }
    }
    #elif EXTAPOLATE_GHOST_CELLS == 2
    // Polynomial extrapolation (chatGPT)
    for (int k = 0; k < nz_ghost; k++)
    {
        for (int j = 0; j < ny; j++)
        {
            array[k][j] = 3*array[nz_ghost][j] - 3*array[nz_ghost+1][j] + array[nz_ghost+2][j];
            array[nz_full-k-1][j] = 3*array[nz_full-nz_ghost-1][j] - 3*array[nz_full-nz_ghost-2][j] + array[nz_full-nz_ghost-3][j];
        }
    }
    #endif
}