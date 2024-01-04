#ifndef FOREGROUND_IO_2D_H__
#define FOREGROUND_IO_2D_H__

#include "data_structures/foreground_data/foreground_data_2D/foreground_variables_struct_2D.h"
#include "data_structures/grid_info/grid_info_2D/grid_info_struct_2D.h"
#include "MPI_module/mpi_info_struct.h"
#include "global_float_precision.h"

void save_foreground_2D(struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, struct MpiInfo *mpi_info, int snap_number, FLOAT_P time);
FLOAT_P load_foreground_2D(struct ForegroundVariables2D *fg, struct GridInfo2D *grid_info, const char *file_path);

#endif // FOREGROUND_IO_2D_H__