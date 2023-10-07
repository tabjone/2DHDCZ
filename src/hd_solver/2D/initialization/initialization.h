#ifndef INITIALIZATION_H
#define INITIALIZATION_H

#include "shared_files.h"
#include "../functions.h"

void initialize_foreground_struct_zeros(struct ForegroundVariables *fg, struct GridInfo *grid_info);
void initialize_foreground_struct_ones(struct ForegroundVariables *fg, struct GridInfo *grid_info);

FLOAT_P gaussian(FLOAT_P x, FLOAT_P y, FLOAT_P x0, FLOAT_P y0, FLOAT_P sigma_x, FLOAT_P sigma_y, FLOAT_P A);

void initialize_foreground_struct_density_pertubation(struct ForegroundVariables *fg, struct BackgroundVariables *bg, struct GridInfo *grid_info);
void initialize_velocity_right(struct ForegroundVariables *fg, struct GridInfo *grid_info);
void initialize_foreground_struct_random(struct ForegroundVariables *fg, struct BackgroundVariables *bg , struct GridInfo *grid_info);

#endif // INITIALIZATION_H