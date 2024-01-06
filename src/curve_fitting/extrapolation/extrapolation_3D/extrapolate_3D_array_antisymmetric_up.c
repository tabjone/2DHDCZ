#include "global_float_precision.h"

void extrapolate_3D_array_antisymmetric_up(FLOAT_P ***array, int nz_full, int nz_ghost, int ny, int nx) 
{
    int start_ghost = nz_full - nz_ghost;

    for (int i = 0; i < nz_ghost; i++) 
    {
        int mirror_row = start_ghost - i - 2;
        for (int j = 0; j < ny; j++) 
        {
            for (int k = 0; k < nx; k++)
            {
                array[start_ghost + i][j][k] = -array[mirror_row][j][k];
            }
        }
    }
}