#ifndef INITIALIZATION_H
#define INITIALIZATION_H

#include "shared_files.h"
#include "../functions.h"
#include "global_parameters.h"

void initialize_foreground(struct ForegroundVariables *fg, struct BackgroundVariables *bg, struct GridInfo *grid_info);

void initialize_foreground_zeros(struct ForegroundVariables *fg, struct GridInfo *grid_info);

FLOAT_P gaussian_2D(FLOAT_P x, FLOAT_P y, FLOAT_P x0, FLOAT_P y0, FLOAT_P sigma_x, FLOAT_P sigma_y, FLOAT_P A);

void initialize_foreground_density_pertubation(struct ForegroundVariables *fg, struct BackgroundVariables *bg, struct GridInfo *grid_info);
void initialize_foreground_entropy_pertubation(struct ForegroundVariables *fg, struct BackgroundVariables *bg, struct GridInfo *grid_info);
void initialize_foreground_velocity_right(struct ForegroundVariables *fg, struct GridInfo *grid_info);
void initialize_foreground_random(struct ForegroundVariables *fg, struct BackgroundVariables *bg , struct GridInfo *grid_info);

void sod_shock_horizontal(struct ForegroundVariables *fg, struct BackgroundVariables *bg, struct GridInfo *grid_info);

#endif // INITIALIZATION_H