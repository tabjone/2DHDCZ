#ifndef FOREGROUND_IO_3D_H__
#define FOREGROUND_IO_3D_H__

#include "data_structures/foreground_data/foreground_data_3D/foreground_variables_struct_3D.h"
#include "data_structures/grid_info/grid_info_3D/grid_info_struct_3D.h"
#include "MPI_module/mpi_info_struct.h"
#include "global_float_precision.h"

void save_foreground_3D(struct ForegroundVariables3D *fg, struct GridInfo3D *grid_info, struct MpiInfo *mpi_info, int snap_number, FLOAT_P time);

#endif // FOREGROUND_IO_3D_H__