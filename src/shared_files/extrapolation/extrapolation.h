#ifndef EXTRAPOLATION_H
#define EXTRAPOLATION_H

#include "shared_files.h"
#include "global_parameters.h"

void extrapolate_2D(struct ForegroundVariables2D *fg, struct GridInfo *grid_info);
void extrapolate_background(struct BackgroundVariables *bg, struct GridInfo *grid_info);
void extrapolate_2D_array(FLOAT_P **array, int nz_full, int nz_ghost, int nx);

#endif // EXTRAPOLATION_H