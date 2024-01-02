#ifndef BOUNDARY_1D_H__
#define BOUNDARY_1D_H__

#include "global_float_precision.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_1D/foreground_variables_struct_1D.h"
#include "data_structures/grid_info/grid_info_1D/grid_info_struct_1D.h"
#include "MPI_module/mpi_info_struct.h"

void apply_vertical_boundary_damping_1D(struct ForegroundVariables1D *fg, struct BackgroundVariables *bg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt);
void calculate_damping_1D(FLOAT_P *damping_factor, struct BackgroundVariables *bg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info);
void update_vertical_boundary_pressure_1D(struct ForegroundVariables1D *fg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info);
void update_vertical_boundary_entropy_velocity_1D(struct ForegroundVariables1D *fg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info);

#endif // BOUNDARY_1D_H__