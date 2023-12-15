#ifndef INITIALIZATION_2D_H
#define INITIALIZATION_2D_H

#include "shared_files.h"
#include "../functions.h"
#include "global_parameters.h"
#include "initialization.h"
#include "global_boundary.h"

#include "../one_time_step/one_time_step.h"
#include "../boundary/boundary.h"

void initialize_foreground(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info);

void initialize_foreground_zeros(struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info);

FLOAT_P gaussian_2D(FLOAT_P x, FLOAT_P y, FLOAT_P x0, FLOAT_P y0, FLOAT_P sigma_x, FLOAT_P sigma_y, FLOAT_P A);

void initialize_foreground_density_pertubation(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info);
void initialize_foreground_entropy_pertubation(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info);
void initialize_foreground_two_entropy_pertubations(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info);
void initialize_foreground_velocity_right(struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info);
void initialize_foreground_random(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg , struct GridInfo2D *grid_info);

void sod_shock_horizontal(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info);
void sod_shock_vertical(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info);

void sod_shock_horizontal_mpi(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info);
void sod_shock_vertical_mpi(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info);

void initialize_foreground_entropy_pertubation_mpi(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info);

void initialization_2D_oscillation_modes(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info);

void new_test_init(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info);

#endif // INITIALIZATION_H