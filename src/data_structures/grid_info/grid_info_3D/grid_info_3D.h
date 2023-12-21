#ifndef GRID_INFO_3D_H
#define GRID_INFO_3D_H

#include "grid_info_struct_3D.h"
#include "../../MPI_module/mpi_info_struct.h"

void allocate_grid_info_struct_3D(struct GridInfo3D **grid_info);
void deallocate_grid_info_struct_3D(struct GridInfo3D *grid_info);
void initialize_grid_info_3D(struct GridInfo3D *grid_info, struct MpiInfo *mpi_info);

#endif // GRID_INFO_3D_H