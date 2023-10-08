#ifndef EXTRAPOLATION_H
#define EXTRAPOLATION_H

#include "shared_files.h"
#include "global_parameters.h"

void extrapolate_2D(struct ForegroundVariables *fg, struct GridInfo *grid_info);
void extrapolate_background(struct BackgroundVariables *bg, struct GridInfo *grid_info);
void extrapolate_2D_array(FLOAT_P **array, int nz_full, int nz_ghost, int ny);


void extrapolate_1D_array_constant_down(FLOAT_P *array, struct GridInfo *grid_info);
void extrapolate_1D_array_constant_up(FLOAT_P *array, struct GridInfo *grid_info);

void extrapolate_2D_array_constant_down(FLOAT_P **array, struct GridInfo *grid_info);
void extrapolate_2D_array_constant_up(FLOAT_P **array, struct GridInfo *grid_info);

void extrapolate_3D_array_constant_down(FLOAT_P ***array, struct GridInfo *grid_info);
void extrapolate_3D_array_constant_up(FLOAT_P ***array, struct GridInfo *grid_info);

#endif // EXTRAPOLATION_H