#ifndef DEALLOCATE_FUNCTIONS_H
#define DEALLOCATE_FUNCTIONS_H

#include "global_float_precision.h"

void deallocate_1D_array(FLOAT_P *array_ptr);
void deallocate_2D_array(FLOAT_P **array_ptr);
void deallocate_3D_array(FLOAT_P ***array_ptr);

#endif // DEALLOCATE_FUNCTIONS_H