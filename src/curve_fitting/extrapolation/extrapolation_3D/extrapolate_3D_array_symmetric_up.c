#include "global_float_precision.h"

void extrapolate_3D_array_symmetric_up(FLOAT_P ***array, int nz_full, int nz_ghost, int ny, int nx) 
{
    for (int i = 0; i < nz_ghost; i++) 
    {
        for (int j = 0; j < ny; j++) 
        {
            for (int k = 0; k < nx; k++)
            {
                array[nz_full - 1 - i][j][k] = array[nz_full - nz_ghost - 1 + i][j][k];
            }
        }
    }
}