#ifndef BOUNDARY_2D_H__
#define BOUNDARY_2D_H__

#include "global_float_precision.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "MPI_module/mpi_info_struct.h"

void apply_vertical_boundary_damping_2D(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt);
void calculate_damping_2D(FLOAT_P *damping_factor, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info);
void update_vertical_boundary_pressure_2D(struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info);
void update_vertical_boundary_entropy_velocity_2D(struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info);

#endif // BOUNDARY_2D_H__