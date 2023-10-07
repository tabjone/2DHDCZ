#include "global_parameters.h"

void copy_3D_array(FLOAT_P ***src, FLOAT_P ***dest, int i_start, int i_end, int j_start, int j_end, int k_start, int k_end)
{
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = j_start; j < j_end; j++)
        {
            for (int k = k_start; k < k_end; k++)
            {
                dest[i][j][k] = src[i][j][k];
            }
        }
    }
}