#ifndef INITIALIZE_FOREGROUND_3D_H_
#define INITIALIZE_FOREGROUND_3D_H_

#include "data_structures/grid_info/grid_info_3D/grid_info_struct_3D.h"
#include "data_structures/foreground_data/foreground_data_3D/foreground_variables_struct_3D.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "MPI_module/MPI_module.h"
#include "global_float_precision.h"

void initialize_foreground_3D(struct ForegroundVariables3D *fg, struct BackgroundVariables *bg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info);

void initialize_foreground_zeros_3D(struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info);


void initialize_foreground_sod_shock_horizontal_3D(struct ForegroundVariables3D *fg, struct BackgroundVariables *bg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info);
void initialize_foreground_sod_shock_vertical_3D(struct ForegroundVariables3D *fg, struct BackgroundVariables *bg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info);

void initialize_foreground_entropy_perturbations_3D(struct ForegroundVariables3D *fg, struct BackgroundVariables *bg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info);

void initialize_foreground_oscillation_modes_3D(struct ForegroundVariables3D *fg, struct BackgroundVariables *bg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info);

void initialize_foreground_random_oscillations_3D(struct ForegroundVariables3D *fg, struct BackgroundVariables *bg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info);

#endif // INITIALIZE_FOREGROUND_3D_H_