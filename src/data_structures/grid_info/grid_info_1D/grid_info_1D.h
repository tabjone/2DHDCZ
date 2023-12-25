#ifndef GRID_INFO_1D_H
#define GRID_INFO_1D_H

#include "grid_info_struct_1D.h"
#include "MPI_module/mpi_info_struct.h"

void allocate_grid_info_struct_1D(struct GridInfo1D **grid_info);
void deallocate_grid_info_struct_1D(struct GridInfo1D *grid_info);
void initialize_grid_info_1D(struct GridInfo1D *grid_info, struct MpiInfo *mpi_info);

#endif // GRID_INFO_2D_H