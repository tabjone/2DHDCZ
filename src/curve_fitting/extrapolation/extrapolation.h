#ifndef EXTRAPOLATION_H
#define EXTRAPOLATION_H

#include "global_float_precision.h"

void extrapolate_1D_array_up(FLOAT_P *array, int nz_full, int nz_ghost);
void extrapolate_1D_array_down(FLOAT_P *array, int nz_ghost);
void extrapolate_1D_array_constant_down(FLOAT_P *array, int nz_ghost);
void extrapolate_1D_array_constant_up(FLOAT_P *array, int nz_full, int nz_ghost);

void extrapolate_2D_array_up(FLOAT_P **array, int nz_full, int nz_ghost, int ny);
void extrapolate_2D_array_down(FLOAT_P **array, int nz_ghost, int ny);
void extrapolate_2D_array_constant_down(FLOAT_P **array, int nz_ghost, int ny);
void extrapolate_2D_array_constant_up(FLOAT_P **array, int nz_full, int nz_ghost, int ny);

void extrapolate_2D_array_antisymmetric_down(FLOAT_P **array, int nz_ghost, int ny);
void extrapolate_2D_array_antisymmetric_up(FLOAT_P **array, int nz_full, int nz_ghost, int ny);

void extrapolate_3D_array_up(FLOAT_P ***array, int nz_full, int nz_ghost, int ny, int nx);
void extrapolate_3D_array_down(FLOAT_P ***array, int nz_ghost, int ny, int nx);
void extrapolate_3D_array_constant_up(FLOAT_P ***array, int nz_full, int nz_ghost, int ny, int nx);
void extrapolate_3D_array_constant_down(FLOAT_P ***array, int nz_ghost, int ny, int nx);

#endif // EXTRAPOLATION_H