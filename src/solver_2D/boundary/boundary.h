#ifndef BOUNDARY_H__
#define BOUNDARY_H__

#include "shared_files.h"
#include "global_parameters.h"

void periodic_boundary_2D(FLOAT_P **array, struct GridInfo2D *grid_info);
void update_vertical_boundary_ghostcells_2D(FLOAT_P **array, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info);
void apply_vertical_boundary_damping(struct ForegroundVariables2D *fg, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt);
void calculate_damping(FLOAT_P *damping_factor, struct BackgroundVariables *bg, struct GridInfo2D *grid_info);
void calculate_damping_mpi(FLOAT_P *damping_factor, struct BackgroundVariables *bg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info);

#endif // BOUNDARY_H__