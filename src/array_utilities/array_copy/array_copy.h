#ifndef ARRAY_COPY_H__
#define ARRAY_COPY_H__

#include "global_float_precision.h"

void copy_1D_array(FLOAT_P *src, FLOAT_P *dest, int i_start, int i_end);
void copy_2D_array(FLOAT_P **src, FLOAT_P **dest, int i_start, int i_end, int j_start, int j_end);
void copy_3D_array(FLOAT_P ***src, FLOAT_P ***dest, int i_start, int i_end, int j_start, int j_end, int k_start, int k_end);

#endif // ARRAY_COPY_H__