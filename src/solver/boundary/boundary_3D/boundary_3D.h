#ifndef BOUNDARY_3D_H__
#define BOUNDARY_3D_H__

#include "global_float_precision.h"
#include "data_structures/background_data/background_variables_struct.h"
#include "data_structures/foreground_data/foreground_data_3D/foreground_variables_struct_3D.h"
#include "data_structures/grid_info/grid_info_3D/grid_info_struct_3D.h"
#include "MPI_module/mpi_info_struct.h"

void apply_vertical_boundary_damping_3D(struct ForegroundVariables3D *fg, struct BackgroundVariables *bg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt);
void calculate_damping_3D(FLOAT_P *damping_factor, struct BackgroundVariables *bg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info);
void update_vertical_boundary_pressure_3D(struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info);
void update_vertical_boundary_entropy_velocity_3D(struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info);

#endif // BOUNDARY_3D_H__