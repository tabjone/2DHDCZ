#ifndef ARRAY_ALLOCATION_ALLOCATE_FUNCTIONS_H
#define ARRAY_ALLOCATION_ALLOCATE_FUNCTIONS_H

void allocate_1D_array(double **array_ptr, int nx);
void allocate_2D_array(double ***array_ptr, int ny, int nx);
void allocate_3D_array(double ****array_ptr, int nz, int ny, int nx);

#endif // ARRAY_ALLOCATION_ALLOCATE_FUNCTIONS_H