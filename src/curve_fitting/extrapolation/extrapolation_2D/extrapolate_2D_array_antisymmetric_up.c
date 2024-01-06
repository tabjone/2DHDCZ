#include "global_float_precision.h"

void extrapolate_2D_array_antisymmetric_up(FLOAT_P **array, int nz_full, int nz_ghost, int ny) 
{

    int start_ghost = nz_full - nz_ghost;

    for (int i = 0; i < nz_ghost; i++) 
    {
        int mirror_row = start_ghost - i - 2;
        for (int j = 0; j < ny; j++) 
        {
            array[start_ghost + i][j] = -array[mirror_row][j];
        }
    }
}