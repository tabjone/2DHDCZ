#include "global_float_precision.h"

void copy_2D_array(FLOAT_P **src, FLOAT_P **dest, int i_start, int i_end, int j_start, int j_end)
{
    /*
    Copies the values from src to dest from i_start to i_end

    Parameters
    ----------
    src : FLOAT_P**
        Pointer to source array
    dest : FLOAT_P**
        Pointer to destination array
    i_start : int
        Index to start copying from in i-direction
    i_end : int
        Index to stop copying at in i-direction
    j_start : int
        Index to start copying from in j-direction
    j_end : int
        Index to stop copying at in j-direction
    */

    for (int i = i_start; i < i_end; i++)
    {
        for (int j = j_start; j < j_end; j++)
        {
            dest[i][j] = src[i][j];
        }
    }
}