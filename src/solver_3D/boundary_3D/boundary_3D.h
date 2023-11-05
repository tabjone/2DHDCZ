#ifndef BOUUNDARY_3D_H__
#define BOUUNDARY_3D_H__

#include "global_parameters.h"
#include "shared_files.h"

void update_vertical_boundary_ghostcells_3D(FLOAT_P ***array, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info);

void periodic_boundary_3D(FLOAT_P ***array, struct GridInfo3D *grid_info);

void calculate_damping_3D(FLOAT_P *damping_factor, struct BackgroundVariables *bg, struct GridInfo3D *grid_info);
void calculate_damping_mpi_3D(FLOAT_P *damping_factor, struct BackgroundVariables *bg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info);

void apply_vertical_boundary_damping_3D(struct ForegroundVariables3D *fg, struct BackgroundVariables *bg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info, FLOAT_P dt);

#endif // BOUUNDARY_3D_H__