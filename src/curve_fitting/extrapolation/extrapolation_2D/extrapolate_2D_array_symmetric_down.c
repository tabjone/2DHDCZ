#include "global_float_precision.h"

void extrapolate_2D_array_symmetric_down(FLOAT_P **array, int nz_ghost, int ny) 
{
    for (int i = 0; i < nz_ghost; i++) 
    {
        for (int j = 0; j < ny; j++) 
        {
            array[i][j] = array[2*nz_ghost-i][j];
        }
    }
}