#ifndef INITIALIZE_FOREGROUND_1D_H
#define INITIALIZE_FOREGROUND_1D_H

#include "data_structures/grid_info/grid_info_1D/grid_info_struct_1D.h"
#include "data_structures/foreground_data/foreground_data_1D/foreground_variables_struct_1D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "MPI_module/MPI_module.h"
#include "global_float_precision.h"

void initialize_foreground_1D(struct ForegroundVariables1D *fg, struct BackgroundVariables *bg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info);

void initialize_foreground_zeros_1D(struct ForegroundVariables1D *fg, struct GridInfo1D *grid_info);


void initialize_foreground_sod_shock_horizontal_1D(struct ForegroundVariables1D *fg, struct BackgroundVariables *bg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info);
void initialize_foreground_sod_shock_vertical_1D(struct ForegroundVariables1D *fg, struct BackgroundVariables *bg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info);

void initialize_foreground_entropy_perturbations_1D(struct ForegroundVariables1D *fg, struct BackgroundVariables *bg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info);

void initialize_foreground_oscillation_modes_1D(struct ForegroundVariables1D *fg, struct BackgroundVariables *bg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info);

void initialize_foreground_random_oscillations_1D(struct ForegroundVariables1D *fg, struct BackgroundVariables *bg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info);

#endif // INITIALIZE_FOREGROUND_1D_H