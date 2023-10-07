#include "global_parameters.h"

void copy_1D_array(FLOAT_P *src, FLOAT_P *dest, int i_start, int i_end)
{
    for (int i = i_start; i < i_end; i++)
    {
        dest[i] = src[i];
    }
}