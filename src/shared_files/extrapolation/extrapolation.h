#ifndef EXTRAPOLATION_H
#define EXTRAPOLATION_H

#include "shared_files.h"

void extrapolate_2D(struct ForegroundVariables2D *fg);
void extrapolate_background(struct BackgroundVariables *bg);
void extrapolate_2D_array(double **array, int nz_full, int nz_ghost, int nx);

#endif // EXTRAPOLATION_H