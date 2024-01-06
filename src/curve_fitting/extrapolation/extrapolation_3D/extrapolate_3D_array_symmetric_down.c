#include "global_float_precision.h"

void extrapolate_3D_array_symmetric_down(FLOAT_P ***array, int nz_ghost, int ny, int nx) 
{
    for (int i = 0; i < nz_ghost; i++) 
    {
        for (int j = 0; j < ny; j++) 
        {
            for (int k = 0; k < nx; k++)
            {
                array[i][j][k] = array[2*nz_ghost-i][j][k];
            }
        }
    }
}