#ifndef EXTRAPOLATION_H
#define EXTRAPOLATION_H

#include "shared_files.h"
#include "global_parameters.h"

void extrapolate_background(struct BackgroundVariables *bg, int nz_full, int nz_ghost, FLOAT_P dz);

void extrapolate_1D_array_up(FLOAT_P *array, int nz_full, int nz_ghost);
void extrapolate_1D_array_down(FLOAT_P *array, int nz_ghost);
void extrapolate_1D_array_constant_down(FLOAT_P *array, int nz_ghost);
void extrapolate_1D_array_constant_up(FLOAT_P *array, int nz_full, int nz_ghost);

void extrapolate_2D_array_up(FLOAT_P **array, int nz_full, int nz_ghost, int ny);
void extrapolate_2D_array_down(FLOAT_P **array, int nz_ghost, int ny);
void extrapolate_2D_array_constant_down(FLOAT_P **array, int nz_ghost, int ny);
void extrapolate_2D_array_constant_up(FLOAT_P **array, int nz_full, int nz_ghost, int ny);

void extrapolate_3D_array_up(FLOAT_P ***array, int nz_full, int nz_ghost, int ny, int nx);
void extrapolate_3D_array_down(FLOAT_P ***array, int nz_ghost, int ny, int nx);
void extrapolate_3D_array_constant_up(FLOAT_P ***array, int nz_full, int nz_ghost, int ny, int nx);
void extrapolate_3D_array_constant_down(FLOAT_P ***array, int nz_ghost, int ny, int nx);

#endif // EXTRAPOLATION_H