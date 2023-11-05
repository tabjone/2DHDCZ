#ifndef FUNCTIONS_3D_H__
#define FUNCTIONS_3D_H__

#include "shared_files.h"
#include "./initialization_3D/initialization_3D.h"
#include "./io_functions_3D/io_functions_3D.h"
#include "./one_time_step_3D/one_time_step_3D.h"

void calculate_grid_info_3D(struct GridInfo3D **grid_info, struct MpiInfo *mpi_info);
void calculate_grid_info_3D_mpi(struct GridInfo3D **grid_info, struct MpiInfo *mpi_info);

#endif // FUNCTIONS_H__