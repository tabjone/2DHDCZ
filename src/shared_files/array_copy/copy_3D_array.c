#include "global_parameters.h"

void copy_3D_array(FLOAT_P ***src, FLOAT_P ***dest, int i_start, int i_end, int j_start, int j_end, int k_start, int k_end)
{
    /*
    Copies the values from src to dest from i_start to i_end

    Parameters
    ----------
    src : FLOAT_P***
        Pointer to source array
    dest : FLOAT_P***
        Pointer to destination array
    i_start : int
        Index to start copying from in i-direction
    i_end : int
        Index to stop copying at in i-direction
    j_start : int
        Index to start copying from in j-direction
    j_end : int
        Index to stop copying at in j-direction
    k_start : int
        Index to start copying from in k-direction
    k_end : int
        Index to stop copying at in k-direction
    */

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