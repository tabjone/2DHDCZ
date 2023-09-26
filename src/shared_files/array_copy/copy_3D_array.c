#include "global_parameters.h"

void copy_3D_array(FLOAT_P ***src, FLOAT_P ***dest, int nz, int ny, int nx)
{
    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nx; k++)
            {
                dest[i][j][k] = src[i][j][k];
            }
        }
    }
}