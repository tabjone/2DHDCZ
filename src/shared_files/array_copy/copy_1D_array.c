#include "global_parameters.h"

void copy_1D_array(FLOAT_P *src, FLOAT_P *dest, int nz)
{
    for (int i = 0; i < nz; i++)
    {
        dest[i] = src[i];
    }
}