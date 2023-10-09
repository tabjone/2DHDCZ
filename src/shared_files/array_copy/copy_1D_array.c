#include "global_parameters.h"

void copy_1D_array(FLOAT_P *src, FLOAT_P *dest, int i_start, int i_end)
{
    /*
    Copies the values from src to dest from i_start to i_end

    Parameters
    ----------
    src : FLOAT_P*
        Pointer to source array
    dest : FLOAT_P*
        Pointer to destination array
    i_start : int
        Index to start copying from in i-direction
    i_end : int
        Index to stop copying at in i-direction
    */

    for (int i = i_start; i < i_end; i++)
    {
        dest[i] = src[i];
    }
}