#ifndef INITIALIZATION_H
#define INITIALIZATION_H

#include "shared_files.h"

void initialize_foreground_struct_zeros(struct ForegroundVariables2D *fg);

double gaussian(double x, double y, double x0, double y0, double sigma_x, double sigma_y, double A);

void initialize_foreground_struct_density_pertubation(struct ForegroundVariables2D *fg);

#endif // INITIALIZATION_H