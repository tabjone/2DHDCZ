#ifndef ARRAY_ALLOCATION_ALLOCATE_FUNCTIONS_H
#define ARRAY_ALLOCATION_ALLOCATE_FUNCTIONS_H

#include "global_parameters.h"


void allocate_1D_array(FLOAT_P **array_ptr, int nx);
void allocate_2D_array(FLOAT_P ***array_ptr, int ny, int nx);
void allocate_3D_array(FLOAT_P ****array_ptr, int nz, int ny, int nx);

#endif // ARRAY_ALLOCATION_ALLOCATE_FUNCTIONS_H