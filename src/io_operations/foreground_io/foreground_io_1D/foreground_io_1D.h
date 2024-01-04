#ifndef FOREGROUND_IO_1D_H__
#define FOREGROUND_IO_1D_H__

#include "data_structures/foreground_data/foreground_data_1D/foreground_variables_struct_1D.h"
#include "data_structures/grid_info/grid_info_1D/grid_info_struct_1D.h"
#include "MPI_module/mpi_info_struct.h"
#include "global_float_precision.h"

void save_foreground_1D(struct ForegroundVariables1D *fg, struct GridInfo1D *grid_info, struct MpiInfo *mpi_info, int snap_number, FLOAT_P time);

#endif // FOREGROUND_IO_1D_H__