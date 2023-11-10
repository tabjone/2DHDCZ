#include "extrapolation.h"

void extrapolate_2D_array_antisymmetric_up(FLOAT_P **array, int nz_full, int nz_ghost, int ny) 
{
    for (int i = 0; i < nz_ghost; i++) 
    {
        for (int j = 0; j < ny; j++) 
        {
            array[nz_full - 1 - i][j] = array[nz_full - nz_ghost - 1 + i][j];
        }
    }
}