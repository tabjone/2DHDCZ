#ifndef INITIALIZE_FOREGROUND_2D_H
#define INITIALIZE_FOREGROUND_2D_H

#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/precalculated_data/precalculated_data_2D/precalculated_data_struct_2D.h"
#include "MPI_module/MPI_module.h"
#include "global_float_precision.h"

void initialize_foreground_2D(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info);

void initialize_foreground_zeros_2D(struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info);


void initialize_foreground_sod_shock_horizontal_2D(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info);
void initialize_foreground_sod_shock_vertical_2D(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info);

void initialize_foreground_entropy_perturbations_2D(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, struct PrecalculatedVariables2D *precalc);

void initialize_foreground_oscillation_modes_2D(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info);

void initialize_foreground_random_oscillations_2D(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info);

#endif // INITIALIZE_FOREGROUND_2D_H